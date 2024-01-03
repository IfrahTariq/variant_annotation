FUNCO_DESC = "Functional annotation from the Funcotator tool."
VEP_CSQ_DESC = "Consequence annotations from Ensembl VEP."
SNPEFF_ANN_DESC = "Functional annotations from SnpEff"
LOF_DESC = "Predicted loss of function effects from VEP"
import pandas as pd
import gzip
import numpy as np
def read_MAF_comments(f):
        description = {}
        colnames = []
        rows = 0
        for l in f:
            l = l.decode("utf-8") if type(l) is not str else l
            if l.startswith("##"):
                rows+=1
                if "FORMAT" in l[:20]:
                    res = l.split("ID=")[1].split(",")[0]
                    desc = l.split("Description=")[1][:-2]
                    description.update({res: desc})
                if "INFO" in l[:20]:
                    res = l.split("ID=")[1].split(",")[0]
                    if res == "FUNCOTATION":
                        print("parsing funcotator special")
                        for val in l.split("Description=")[1][:-2].split("|"):
                            val = val.split("Funcotation fields are: ")[-1]
                            description.update({val: FUNCO_DESC})
                    elif res == "ANN":
                        print("parsing funcitonal annotation from SnpEff")
                        l = l.replace(" / ", "_").replace(".", "_")
                        for val in l.split("Description=")[1][:-5].split(" | "):
                            val = "snpeff_" + val.split("Functional annotations: '")[-1]
                            description.update({val: SNPEFF_ANN_DESC})
                    elif res == "LOF":
                        print("parsing predicted LOF status from SnpEff")
                        for val in l.split("Description=")[1][:-4].split(" | "):
                            val = "lof_" + val.split("Format: '")[-1]
                            description.update({val: LOF_DESC})
                    elif res == "CSQ":
                        print("parsing VEP CSQ")
                        for val in l.split("Description=")[1][:-3].split("|"):
                            val = "vep_" + val.split("Format: ")[-1]
                            description.update({val: VEP_CSQ_DESC})
                    elif res == "REF":
                        description.update(
                            {"REF_FLAG": l.split("Description=")[1][:-2]}
                        )
                    else:
                        desc = l.split("Description=")[1][:-2]
                        description.update({res: desc})
            elif l.startswith("Hugo_Symbol"):
                rows+=1
                colnames = l.split("\t")
        return description, colnames, rows

def read_VCF_comments(f):
    description = {}
    colnames = []
    rows = 0
    filters = []
    formats = []
    infos = []
    for l in f:
        l = l.decode("utf-8") if type(l) is not str else l
        if l.startswith("##"):
            rows+=1
            if "FILTER=" in l:
                res = l.split("ID=")[1].split(",")[0]
                desc = l.split("Description=")[1][:-2]
                description.update({res: desc})
                filters.append(res)
            if "FORMAT" in l:
                res = l.split("ID=")[1].split(",")[0]
                desc = l.split("Description=")[1][:-2]
                description.update({res: desc})
                formats.append(res)
            if "INFO" in l:
                res = l.split("ID=")[1].split(",")[0]
                # print(res)
                if res == "FUNCOTATION":
                    print("parsing funcotator special")
                    for val in l.split("Description=")[1][:-2].split("|"):
                        val = val.split("Funcotation fields are: ")[-1]
                        description.update({val: FUNCO_DESC})
                elif res == "ANN":
                    print("parsing funcitonal annotation from SnpEff")
                    l = l.replace(" / ", "_").replace(".", "_")
                    for val in l.split("Description=")[1][:-5].split(" | "):
                        val = "snpeff_" + val.split("Functional annotations: '")[-1]
                        description.update({val: SNPEFF_ANN_DESC})
                elif res == "LOF":
                    print("parsing predicted LOF status from SnpEff")
                    for val in l.split("Description=")[1][:-4].split(" | "):
                        val = "lof_" + val.split("Format: '")[-1]
                        description.update({val: LOF_DESC})
                elif res == "CSQ":
                    print("parsing VEP CSQ")
                    for val in l.split("Description=")[1][:-3].split("|"):
                        val = "vep_" + val.split("Format: ")[-1]
                        description.update({val: VEP_CSQ_DESC})
                elif res == "REF":
                    description.update(
                        {"REF_FLAG": l.split("Description=")[1][:-2]}
                    )
                infos.append(res)
        elif l.startswith("#CHROM"):
            colnames = l[1:-1].split("\t")
        else:
            break
    return description, colnames, rows

def funcotate_and_vep(funcotator_maf_file, vep_vcf_file):
    f = open(funcotator_maf_file, "r")
    description, colnames, rows = read_MAF_comments(f)
    maf = pd.read_csv(funcotator_maf_file, sep="\t", skiprows=rows)

    f1 = gzip.open(vep_vcf_file, "rb")
    description, colnames, rows = read_VCF_comments(f1)
    vcf = pd.read_csv(vep_vcf_file, compression='gzip', sep="\t", skiprows=rows)

    vcf_main = vcf[['#CHROM', 'POS', 'ID','REF','ALT','QUAL','FILTER']]
    # Split the INFO column by ";"
    info_list = vcf['INFO'].str.split(';')
    # Split each key-value pair by "=" and create a new DataFrame
    vcf_info = pd.DataFrame([dict(item.split('=') for item in sublist if item.count('=') == 1) for sublist in info_list])

    csq_df = vcf_info['CSQ'].str.split('|', expand=True)
    csq_cols = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|pLI_gene_value|LoFtool".split("|")
    try:
        csq_df.columns = csq_cols
    except:
        csq_df.columns = [f'ANN_fix{x}' for x in range(len(csq_df.columns.tolist()))]

    ann_df = vcf_info['ANN'].str.split('|', expand=True)
    ann_cols = 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'.split("|")
    # Drop columns where the majority of values are None
    ann_df = ann_df.dropna(thresh= len(ann_df) // 2, axis=1)
    try:
        ann_df.columns = ann_cols
    except:
        ann_df.columns = [f'ANN_fix{x}' for x in range(len(ann_df.columns.tolist()))]

    vcf_info = vcf_info.drop(['ANN', 'CSQ'], axis=1)

    vcf_normal = vcf.iloc[:,-1].str.replace(",", "-").str.split(":")
    vcf_normal = pd.DataFrame(vcf_normal.tolist()).iloc[:,:8]
    try:
        vcf_normal.columns = [f'NORMAL_{x}' for x in vcf.FORMAT[0].split(":")]
    except:
        vcf_normal.columns = [f'NORMAL_fix{x}' for x in range(len(vcf_normal.columns.tolist()))]
    vcf_tumor = vcf.iloc[:,-2].str.replace(",", "-").str.split(":")
    vcf_tumor = pd.DataFrame(vcf_tumor.tolist()).iloc[:,:8]
    try:
        vcf_tumor.columns = [f'TUMOR_{x}' for x in vcf.FORMAT[0].split(":")]
    except:
        vcf_tumor.columns = [f'TUMOR_fix{x}' for x in range(len(vcf_tumor.columns.tolist()))]
        
    vcf_final = pd.concat([vcf_main, vcf_info, csq_df, ann_df, vcf_tumor, vcf_normal], axis=1)
    
    maf = maf.set_index('Start_Position')
    vcf_final = vcf_final.set_index('POS')

    common_cols = list(set(maf.columns.tolist()).intersection(vcf_final.columns.tolist()))
    maf_dif_cols = list(set(maf.columns.tolist()).difference(vcf_final.columns.tolist()))
    vcf_dif_cols = list(set(vcf_final.columns.tolist()).difference(maf.columns.tolist()))
    maf_unique = maf[maf_dif_cols]
    vcf_unique = vcf_final[common_cols + vcf_dif_cols]
    final_df = pd.merge(maf_unique, vcf_unique, left_index=True, right_index=True)
    final_df = final_df.dropna(how='all', axis=1)

    return final_df, maf, vcf_final

def cleanup_vcf(vcf):
    replace_empty=REPLACE_EMPTY
    print("replacing empty characters:")
    vcf = vcf.replace(replace_empty)
    ## determine which variants are associated with splicing
    try:
        subvcf = vcf[(vcf["OC_spliceai__ds_ag"] != "")]
        loc = subvcf[
            (subvcf["OC_spliceai__ds_ag"].astype(float) >= 0.5)
            | (subvcf["OC_spliceai__ds_al"].astype(float) >= 0.5)
            | (subvcf["OC_spliceai__ds_dg"].astype(float) >= 0.5)
            | (subvcf["OC_spliceai__ds_dl"].astype(float) >= 0.5)
        ].index
        vcf['assoc_splicing'] = 0
        vcf.loc[loc, "assoc_splicing"] = 1
    except:
        vcf['assoc_splicing'] = 'fix'
    ## remain OpenCravat columns
    vcf = vcf.rename(columns=RENAME_OC)
    ## drop gnomad-genome columns
    vcf = vcf.drop(columns = vcf[[col for col in vcf.columns if 'gnomADg' in col]].columns.tolist())
    
    if "PolyPhen" in vcf.columns.tolist():
        ## cleanup PolyPhen and SIFT columns
        vcf_phen = vcf[[col for col in vcf.columns if 'PolyPhen' in col]]
        # Get the part before "(" if it exists
        vcf["PolyPhen_anno"] = [x.split("(")[0] if "(" in x else x for x in vcf_phen.PolyPhen.tolist()]
        # Get the part after "(" if it exists
        numbers = [x.split("(")[1] if "(" in x else '' for x in vcf_phen.PolyPhen.tolist()]
        vcf["PolyPhen_score"] = [x.replace(')', '') for x in numbers]
    if "SIFT" in vcf.columns.tolist():
        vcf_sift = vcf[[col for col in vcf.columns if 'SIFT' in col]]#.columns.tolist()
        # Get the part before "(" if it exists
        vcf["SIFT_anno"] = [x.split("(")[0] if "(" in x else x for x in vcf_sift.SIFT.tolist()]
        # Get the part after "(" if it exists
        numbers = [x.split("(")[1] if "(" in x else '' for x in vcf_sift.SIFT.tolist()]
        vcf["SIFT_score"] = [x.replace(')', '') for x in numbers]
    return vcf


RENAME_OC = {
 'OC_base__achange': 'base_achange',
 'OC_base__coding': 'base_coding',
 'OC_base__so': 'base_so',
 'OC_base__all_mappings': 'base_all_mappings',
 'OC_base__cchange': 'base_cchange',
 'OC_base__transcript': 'base_transcript',
 'OC_base__hugo': 'base_hugo',
 'OC_alfa__total_freq': 'alfa_total_freq',
 'OC_alfa__total_alt': 'alfa_total_alt',
 'OC_ccre_screen___group': 'encode_group',
 'OC_ccre_screen__bound': 'encode_bound',
 'OC_ccre_screen__acc_d': 'encode_acc_d',
 'OC_ccre_screen__acc_e': 'encode_acc_e',
 'OC_spliceai__ds_dl': 'spliceai_ds_dl',
 'OC_spliceai__dp_al': 'spliceai_dp_al',
 'OC_spliceai__ds_ag': 'spliceai_ds_ag',
 'OC_spliceai__ds_al': 'spliceai_ds_al',
 'OC_spliceai__dp_dg': 'spliceai_dp_dg',
 'OC_spliceai__ds_dg': 'spliceai_ds_dg',
 'OC_spliceai__dp_dl': 'spliceai_ds_dl',
 'OC_spliceai__dp_ag': 'spliceai_dp_ag',
 'OC_pharmgkb__id': 'pharmgkb_id',
 'OC_pharmgkb__chemicals': 'pharmgkb_chemicals',
 'OC_pharmgkb__drug_assoc': 'pharmgkb_drug_assoc',
 'OC_pharmgkb__pheno_cat': 'pharmgkb_pheno_cat',
 'OC_gwas_catalog__or_beta': 'gwas_or_beta',
 'OC_gwas_catalog__disease': 'gwas_disease',
 'OC_gwas_catalog__ci': 'gwas_ci',
 'OC_gwas_catalog__pmid': 'gwas_pmid',
 'OC_gwas_catalog__init_samp': 'gwas_init_samp',
 'OC_gwas_catalog__rep_samp': 'gwas_rep_samp',
 'OC_gwas_catalog__pval': 'gwas_pval',
 'OC_gwas_catalog__risk_allele': 'gwas_risk_allele',
 'OC_provean__score': 'provean_score',
 'OC_provean__rankscore': 'provean_rankscore',
 'OC_provean__all': 'provean_all',
 'OC_provean__transcript': 'provean_transcript',
 'OC_provean__uniprot': 'provean_uniprot',
 'OC_provean__prediction': 'provean_prediction',
 'OC_revel__all': 'revel_all',
 'OC_revel__score': 'revel_score',
 'OC_revel__transcript': 'revel_transcript',
 'OC_revel__rankscore': 'revel_rankscore',
 'OC_original_input__ref_base': 'original_input_ref_base',
 'OC_original_input__chrom': 'original_input_chrom',
 'OC_original_input__alt_base': 'original_input_alt',
 'OC_original_input__pos': 'original_input_pos',
 'OC_gtex__gtex_gene': 'gtex_gene',
 'OC_gtex__gtex_tissue': 'gtex_tissue',
 'OC_hess_drivers__is_driver': 'hess_driver',
 'OC_hess_drivers__signature': 'hess_signture',
 }

REPLACE_EMPTY = {
    np.nan: "",
    None: "",
    "Unknown": "",
    "unknown": "",
    "UNKNOWN": "",
    ",": "",
    ",,": "",
    ",,,": "",
    ",%3B,": "",
    ",;,": "",
    ",%2C,": "",
    ",,,,": "",
    ",,,,,": "",
    ",,,,,,": "",
    ",;,,;,": "",
    ",%3B,,%3B,": "",
    ",%2C,,%2C,": "",
    ",,,,,,,": "",
    ",,,,,,,,": "",
    ",,,,,,,,,": "",
    ",;,,;,,;,": "",
    ",%3B,,%3B,,%3B,": "",
    ",%2C,,%2C,,%2C,": "",
    "|": "",
    "||": "",
    "|||": "",
    "||||": "",
    "|||||": "",
    "||||||": "",
    "Approved": "Y",
    "yes": "Y",
    "Yes": "Y",
    "YES": "Y",
    "y": "Y",
    "True": "Y",
    "False": "",
    True: "Y",
    False: "",
    "true": "Y",
    "false": "",
}
