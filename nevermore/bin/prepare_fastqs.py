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
	pairs = {}
	for p1, p2 in itertools.product(
		tuple(p[:-1] for p, _ in r1),
		tuple(p[:-1] for p, _ in r2)
	):
		pairs.setdefault(p1, set()).add(1)
		pairs.setdefault(p2, set()).add(2)

	for p, counts in pairs.items():
		if len(counts) < 2:
			raise ValueError(f"Missing mates for prefix {p}, files={str(counts)}")
		elif len(counts) > 2:
			raise ValueError(f"Too many files for prefix {p}, files={str(counts)}")
		else:
			...


def transfer_file(source, dest, remote_input=False):
	if not source.name.endswith(".gz"):
		with open(dest, "wt") as _out:
			subprocess.run(("gzip", "-c", source.resolve()), stdout=_out)
	elif remote_input:
		shutil.copyfile(source.resolve(), dest)
	else:
		pathlib.Path(dest).symlink_to(source.resolve())


def transfer_multifiles(files, dest, remote_input=False, gzipped=False):
	if len(files) > 1:
		src_files = tuple(f.resolve() for f in files)
		cat_cmd = ("cat", ) + src_files
		if not gzipped:
			cat_pr = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
			with open(dest, "wt") as _out:
				subprocess.run(("gzip", "-c", "-"), stdin=cat_pr.stdout, stdout=_out)
		else:
			with open(dest, "wt") as _out:
				subprocess.run(cat_cmd, stdout=_out)
	elif not gzipped:
		with open(dest, "wt") as _out:
			subprocess.run(("gzip", "-c",) + src_files, stdout=_out)
	elif remote_input:
		shutil.copyfile(files[0].resolve(), dest)
	else:
		pathlib.Path(dest).symlink_to(files[0].resolve())


def process_sample(sample, fastqs, output_dir, remove_suffix=None, remote_input=False):
	
	if len(fastqs) == 1:
		sample_sub = re.sub(r"[._]singles?", "", sample)
		if sample_sub != sample:
			sample = sample_sub + ".singles"
		sample_dir = os.path.join(output_dir, sample)
		pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)

		dest = os.path.join(sample_dir, f"{sample}_R1.fastq.gz")
		transfer_file(fastqs[0], dest, remote_input=remote_input)

	else:

		gzips = {f for f in fastqs if f.name.endswith(".gz")}
		no_gzips = {f for f in fastqs if not f.name.endswith(".gz")}
		if gzips and no_gzips:
			raise ValueError(f"sample: {sample} has mixed gz and uncompressed input files. Please check.")

		prefixes = [re.sub(r"\.(fastq|fq).gz$", "", os.path.basename(f.name)) for f in fastqs]
		if remove_suffix:
			prefixes = [re.sub(remove_suffix + r"$", "", p) for p in prefixes]

		print("PRE", prefixes, file=sys.stderr)

		r1 = [(p, f) for p, f in zip(prefixes, fastqs) if p.endswith("1")]
		r2 = [(p, f) for p, f in zip(prefixes, fastqs) if p.endswith("2")]
		others = set(fastqs).difference({f for _, f in r1}).difference({f for _, f in r2})

		assert len(r2) == 0 or len(r1) == len(r2), "R1/R2 sets are not of the same length"
		check_pairwise(r1, r2)

		r1 = sorted(f for _, f in r1)
		r2 = sorted(f for _, f in r2)

		print("R1", r1, file=sys.stderr)
		print("R2", r2, file=sys.stderr)

		sample_dir = os.path.join(output_dir, sample)
		pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)

		if r1:
			dest = os.path.join(sample_dir, f"{sample}_R1.fastq.gz")
			transfer_multifiles(r1, dest, remote_input=remote_input, gzipped=bool(gzips))
		if r2:
			dest = os.path.join(sample_dir, f"{sample}_R2.fastq.gz")
			transfer_multifiles(r2, dest, remote_input=remote_input, gzipped=bool(gzips))
		if others:
			sample_dir = sample_dir + ".singles"
			pathlib.Path(sample_dir).mkdir(parents=True, exist_ok=True)
			dest = os.path.join(sample_dir, f"{sample}.singles_R1.fastq.gz")
			transfer_multifiles(others, dest, remote_input=remote_input, gzipped=bool(gzips))
		

def is_fastq(f):
	prefix, suffix = os.path.splitext(f)
	if suffix in ("fastq", "fq"):
		return True
	if suffix == ".gz":
		_, suffix = os.path.splitext(prefix)
		return suffix in ("fastq", "fq")
	return False


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input_dir", type=str, default=".")
	ap.add_argument("-o", "--output_dir", type=str, default="prepared_samples")
	ap.add_argument("--remote-input", action="store_true")
	ap.add_argument("--remove-suffix", type=str, default=None)

	args = ap.parse_args()

	pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)
	
	fastqs = sorted(
		f 
		for f in os.listdir(args.input_dir)
		if is_fastq(f)
	)
	assert fastqs, f"Could not find any fastq files in {args.input_dir}."

	samples = {}

	for f in fastqs:
		full_f = pathlib.Path(os.path.join(args.input_dir, f))
		if full_f.is_symlink():
			link_target = full_f.resolve()
			sample = os.path.basename(os.path.dirname(link_target))
			if not sample:
				raise NotImplementedError("Flat-directories not implemented.")
			samples.setdefault(sample, []).append(full_f)

	for sample, fastqs in samples.items():
		try:
			process_sample(
				sample, fastqs, args.output_dir,
				remove_suffix=args.remove_suffix, remote_input=args.remote_input
			)
		except Exception as e:
			raise ValueError(f"Encountered problems processing sample '{args.sample_id}': {e}.\nPlease check your file names.")


if __name__ == "__main__":
	main()