import utils
import pandas as pd
import numpy as np
import argparse


def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("vep_file")
    parser.add_argument("funcotator_file")
    parser.add_argument("sample_name", nargs="?", default=None)
    args = parser.parse_args()

    vep_file = args.vep_file
    funcotator_file = args.funcotator_file
    sample_name = args.sample_name

    print(
        "inputs: vcf_file:",
        vep_file,
        ", funcotator_file:",
        funcotator_file,
        ", sample_name:",
        sample_name
    )
    vcf, _, _ = utils.funcotate_and_vep(funcotator_file, vep_file)
    vcf = utils.cleanup_vcf(vcf)
    vcf.to_csv(sample_name + "_vep_funcotator.csv")
    print(f"finished, saved csv")

if __name__ == "__main__":
    main()
