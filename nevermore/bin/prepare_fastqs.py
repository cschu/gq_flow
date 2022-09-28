#!/usr/bin/env python3

import argparse
import itertools
import os
import pathlib
import re
import shutil
import subprocess
import sys


def check_pairwise(r1, r2):
	""" Checks if two sets of read files contain the same prefixes.

	Input:
	 - r1/r2: lists of tuples (prefix, filename) of R1/R2 paired-end read files

	Raises error if one prefix is not found in either r1 or r2.

	"""
	r1_r2 = tuple(prefix[:-1] for prefix, _ in itertools.chain(r1, r2))
	for prefix in set(r1_r2):
		if r1_r2.count(prefix) != 2:
			raise ValueError(f"Missing mates for prefix {prefix}.")


def transfer_file(source, dest, remote_input=False):
	""" Transfers a file depending on its location and compressed state.

	Input:
	 - path to source file
	 - path to destination file
	 - whether source file is considered to be located on a remote file system

	"""
	if not source.name.endswith(".gz"):
		# if file is not gzipped, gzip it to destination
		with open(dest, "wt") as _out:
			subprocess.run(("gzip", "-c", source.resolve()), stdout=_out)
	elif remote_input:
		# if file is on remote file system, copy it to destination
		shutil.copyfile(source.resolve(), dest)
	else:
		# if file is gzipped and on local fs, just symlink it
		pathlib.Path(dest).symlink_to(source.resolve())


def transfer_multifiles(files, dest, remote_input=False, gzipped=False):
	""" Transfers a set of files depending on their location and compressed state.

	Input:
	 - list of source file paths
	 - path to destination file
	 - whether source files are considered to be located on a remote file system
	 - whether source files are gzipped

	"""
	if len(files) > 1:
		src_files = tuple(f.resolve() for f in files)
		cat_cmd = ("cat", ) + src_files
		if not gzipped:
			# multiple uncompressed files will be cat | gzipped
			cat_pr = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
			with open(dest, "wt") as _out:
				subprocess.run(("gzip", "-c", "-"), stdin=cat_pr.stdout, stdout=_out)
		else:
			# multiple compressed files with just be concatenated
			with open(dest, "wt") as _out:
				subprocess.run(cat_cmd, stdout=_out)
	else:
		transfer_file(files[0], dest, remote_input=remote_input)


def process_sample(sample, fastqs, output_dir, remove_suffix=None, remote_input=False):
	""" Checks if a set of fastq files in a directory is a valid collection
	and transfers files to a destination dir upon success.

	Input:
	 - sample_id
	 - list of fastq files
	 - path to output directory
	 - suffix to strip off from filenames (e.g. _001)
	 - fastq files are located on remote file system
	"""

	if len(fastqs) == 1:
		# remove potential "single(s)" string from single fastq file name prefix
		sample_sub = re.sub(r"[._]singles?", "", sample)
		if sample_sub != sample:
			# if file name had single(s) pattern,
			# attach it to the end of the prefix
			sample = sample_sub + ".singles"
		sample_dir = os.path.join(output_dir, sample)
		pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)

		dest = os.path.join(sample_dir, f"{sample}_R1.fastq.gz")
		transfer_file(fastqs[0], dest, remote_input=remote_input)

	else:

		# check if all fastq files are either gzipped or not
		gzips = {f for f in fastqs if f.name.endswith(".gz")}
		no_gzips = {f for f in fastqs if not f.name.endswith(".gz")}
		if gzips and no_gzips:
			raise ValueError(f"sample: {sample} has mixed gz and uncompressed input files. Please check.")

		# extract the file name prefixes
		prefixes = [re.sub(r"\.(fastq|fq).gz$", "", os.path.basename(f.name)) for f in fastqs]
		if remove_suffix:
			# remove suffix pattern if requested
			prefixes = [re.sub(remove_suffix + r"$", "", p) for p in prefixes]

		print("PRE", prefixes, file=sys.stderr)

		# partition fastqs into R1, R2, and 'other' sets
		r1 = [(p, f) for p, f in zip(prefixes, fastqs) if re.search(r"[._R]1$", p)]
		r2 = [(p, f) for p, f in zip(prefixes, fastqs) if re.search(r"[._R]2$", p)]
		others = list(set(fastqs).difference({f for _, f in r1}).difference({f for _, f in r2}))

		# check if R1/R2 sets have equal sizes or are empty
		# R1 empty: potential scRNAseq (or any protocol with barcode reads in R1)
		# R2 empty: typical single end reads with (R?)1 suffix
		assert len(r2) == 0 or len(r1) == 0 or (r1 and len(r1) == len(r2)), "R1/R2 sets are not of the same length"

		# if R1 and R2 are of equal size, check if the prefixes match
		if len(r1) == len(r2) and r1:
			check_pairwise(r1, r2)

		# sort R1/R2 for concatenation, get rid off prefixes
		r1 = sorted(f for _, f in r1)
		r2 = sorted(f for _, f in r2)

		print("R1", r1, file=sys.stderr)
		print("R2", r2, file=sys.stderr)

		sample_dir = os.path.join(output_dir, sample)
		pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)

		if r1:
			# if R1 is not empty, transfer R1-files
			dest = os.path.join(sample_dir, f"{sample}_R1.fastq.gz")
			transfer_multifiles(r1, dest, remote_input=remote_input, gzipped=bool(gzips))
		if r2:
			# if R2 is not empty, transfer R2-files,
			# if R1 is empty, rename R2 to R1 so that files can be processed as normal single-end
			target_r = "R2" if r1 else "R1"
			dest = os.path.join(sample_dir, f"{sample}_R2.fastq.gz")
			transfer_multifiles(r2, dest, remote_input=remote_input, gzipped=bool(gzips))
		if others:
			# if single-end reads exist,
			# transfer them to <sample>.singles
			# these will be processed independently and merged with the paired-end reads
			# at a later stage
			sample_dir = sample_dir + ".singles"
			pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)
			dest = os.path.join(sample_dir, f"{sample}.singles_R1.fastq.gz")
			transfer_multifiles(others, dest, remote_input=remote_input, gzipped=bool(gzips))
		

def is_fastq(f):
	""" Checks if a file is a fastq file (compressed or uncompressed.)

	Input:
	 - filename

	Output:
	 - true if file is fastq else false

	"""
	prefix, suffix = os.path.splitext(f)
	if suffix in (".fastq", ".fq"):
		return True
	if suffix == ".gz":
		_, suffix = os.path.splitext(prefix)
		return suffix in (".fastq", ".fq")
	return False


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input_dir", type=str, default=".")
	ap.add_argument("-o", "--output_dir", type=str, default="prepared_samples")
	ap.add_argument("-p", "--prefix", type=str, required=True)
	ap.add_argument("--remote-input", action="store_true")
	ap.add_argument("--remove-suffix", type=str, default=None)

	args = ap.parse_args()

	pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)
	
	# collect all fastq files from input directory
	# assumption: fastq files are sym-linked into input_dir from a sample/files directory tree
	# i.e., from a nextflow.Channel()
	fastqs = sorted(
		f 
		for f in os.listdir(args.input_dir)
		if is_fastq(f)
	)
	assert fastqs, f"Could not find any fastq files in '{args.input_dir}'."

	samples = {}

	# resolve the symlinks and group by sample ids
	for f in fastqs:
		full_f = pathlib.Path(os.path.join(args.input_dir, f))
		if full_f.is_symlink():
			link_target = full_f.resolve()
			sample, *fpath = str(link_target).replace(args.prefix, "").lstrip("/").split("/")

			if not fpath:
				raise NotImplementedError("Flat-directories not implemented.")
			samples.setdefault(sample, []).append(full_f)

	# check and transfer the files
	for sample, fastqs in samples.items():
		try:
			process_sample(
				sample, fastqs, args.output_dir,
				remove_suffix=args.remove_suffix, remote_input=args.remote_input
			)
		except Exception as e:
			raise ValueError(f"Encountered problems processing sample '{sample}': {e}.\nPlease check your file names.")


if __name__ == "__main__":
	main()