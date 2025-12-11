#!/usr/bin/env python3
"""
Create CSV files for each viral family with contig cohort information
Based on the R script plot_contig_heatmaps.R
"""

import pandas as pd
from pathlib import Path

# Define the base directory
base_dir = Path("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs")

# Define all contigs and their cohort presence based on the R script
viral_families = {
    'Astroviridae': {
        'all_contigs': ["contig_3357", "contig_34440", "contig_3567", "contig_3669",
                       "contig_43336", "contig_4647", "contig_49047", "contig_49065", "contig_8153"],
        'US_contigs': ["contig_34440", "contig_3567", "contig_4647", "contig_49065"],
        'Italy_contigs': ["contig_34440", "contig_3567", "contig_4647", "contig_49065"]
    },

    'Anelloviridae': {
        'all_contigs': ["contig_30686", "contig_14733", "contig_42538", "contig_7638", "contig_39636",
                       "contig_8637", "contig_10228", "contig_9781", "contig_7962", "contig_23715",
                       "contig_36641", "contig_13000", "contig_41090", "contig_25647", "contig_6273",
                       "contig_30674", "contig_34347", "contig_21558", "contig_6488", "contig_33822",
                       "contig_10424", "contig_41802", "contig_7169", "contig_22197", "contig_25919",
                       "contig_8612", "contig_21645", "contig_26562", "contig_46055", "contig_15372",
                       "contig_38908", "contig_47848", "contig_49946", "contig_14012", "contig_22018",
                       "contig_9359", "contig_25612", "contig_11020", "contig_7799", "contig_33967",
                       "contig_49814", "contig_40383", "contig_40400", "contig_49639", "contig_8416",
                       "contig_20459", "contig_11215", "contig_35851", "contig_36742", "contig_45997",
                       "contig_13393", "contig_14584", "contig_8872", "contig_9163", "contig_9685",
                       "contig_46929", "contig_26624", "contig_43766", "contig_9281", "contig_11113",
                       "contig_9335", "contig_5833", "contig_685", "contig_28092", "contig_26972",
                       "contig_10357", "contig_27323", "contig_8313", "contig_36473", "contig_30210",
                       "contig_17146", "contig_23561", "contig_7355", "contig_8111", "contig_8203",
                       "contig_35169", "contig_16851", "contig_6571", "contig_22898", "contig_10173",
                       "contig_32356", "contig_14117", "contig_6711", "contig_17385", "contig_40927",
                       "contig_7778", "contig_25715", "contig_29277", "contig_35520", "contig_38987",
                       "contig_9065", "contig_11182", "contig_17660", "contig_18632", "contig_7587",
                       "contig_42090", "contig_46346", "contig_33693", "contig_8484", "contig_8505",
                       "contig_12655", "contig_30042", "contig_7798", "contig_12607", "contig_7294",
                       "contig_42211", "contig_11040", "contig_13635", "contig_6399", "contig_7262",
                       "contig_8253", "contig_12364", "contig_686", "contig_8312", "contig_6829",
                       "contig_30191", "contig_20764", "contig_38226", "contig_17021", "contig_10247",
                       "contig_10528", "contig_14516", "contig_35806", "contig_13875", "contig_44065",
                       "contig_7525", "contig_7210", "contig_40324", "contig_45257", "contig_25052",
                       "contig_14022", "contig_19233", "contig_33077", "contig_7521", "contig_7219",
                       "contig_11970", "contig_32845", "contig_11420", "contig_7949", "contig_11066",
                       "contig_35170", "contig_24801", "contig_13532", "contig_9099", "contig_6926",
                       "contig_15297", "contig_7325", "contig_7956", "contig_9565", "contig_8930", "contig_12325"],
        'US_contigs': ["contig_42538", "contig_7638", "contig_39636", "contig_8637", "contig_9781",
                      "contig_23715", "contig_36641", "contig_6273", "contig_34347", "contig_21558",
                      "contig_33822", "contig_22197", "contig_25919", "contig_26562", "contig_46055",
                      "contig_49946", "contig_14012", "contig_25612", "contig_11020", "contig_7799",
                      "contig_33967", "contig_49814", "contig_40383", "contig_40400", "contig_49639",
                      "contig_8416", "contig_11215", "contig_35851", "contig_36742", "contig_45997",
                      "contig_13393", "contig_14584", "contig_46929", "contig_26624", "contig_43766",
                      "contig_11113", "contig_9335", "contig_5833", "contig_685", "contig_28092",
                      "contig_10357", "contig_36473", "contig_17146", "contig_8203", "contig_35169",
                      "contig_6571", "contig_22898", "contig_10173", "contig_6711", "contig_7778",
                      "contig_25715", "contig_29277", "contig_38987", "contig_9065", "contig_11182",
                      "contig_18632", "contig_7587", "contig_42090", "contig_46346", "contig_33693",
                      "contig_8484", "contig_7798", "contig_7294", "contig_42211", "contig_13635",
                      "contig_8253", "contig_12364", "contig_686", "contig_8312", "contig_6829",
                      "contig_38226", "contig_17021", "contig_10247", "contig_14516", "contig_13875",
                      "contig_44065", "contig_7525", "contig_7210", "contig_45257", "contig_25052",
                      "contig_19233", "contig_33077", "contig_7521", "contig_11970", "contig_32845",
                      "contig_7949", "contig_35170", "contig_24801", "contig_13532", "contig_7325",
                      "contig_7956"],
        'Italy_contigs': ["contig_30686", "contig_14733", "contig_9781", "contig_7962", "contig_23715",
                         "contig_13000", "contig_41090", "contig_25647", "contig_30674", "contig_6488",
                         "contig_10424", "contig_41802", "contig_7169", "contig_8612", "contig_21645",
                         "contig_46055", "contig_15372", "contig_38908", "contig_47848", "contig_49946",
                         "contig_9359", "contig_49814", "contig_20459", "contig_11215", "contig_45997",
                         "contig_13393", "contig_8872", "contig_9163", "contig_9685", "contig_9281",
                         "contig_11113", "contig_685", "contig_26972", "contig_27323", "contig_8313",
                         "contig_30210", "contig_23561", "contig_8111", "contig_16851", "contig_32356",
                         "contig_14117", "contig_17385", "contig_40927", "contig_35520", "contig_17660",
                         "contig_46346", "contig_33693", "contig_8484", "contig_8505", "contig_12655",
                         "contig_30042", "contig_7798", "contig_12607", "contig_11040", "contig_6399",
                         "contig_686", "contig_30191", "contig_20764", "contig_10528", "contig_35806",
                         "contig_13875", "contig_44065", "contig_40324", "contig_14022", "contig_7521",
                         "contig_7219", "contig_11420", "contig_11066", "contig_9099", "contig_6926",
                         "contig_15297", "contig_9565", "contig_8930"]
    },

    'Adenoviridae': {
        'all_contigs': ["contig_48716", "contig_2664", "contig_2451", "contig_2453", "contig_48715",
                       "contig_13045", "contig_2450", "contig_2464", "contig_33906", "contig_2452",
                       "contig_5040", "contig_8600", "contig_7385", "contig_48956"],
        'US_contigs': ["contig_2664", "contig_2464", "contig_33906"],
        'Italy_contigs': ["contig_48716", "contig_2664", "contig_2451", "contig_2453", "contig_48715",
                         "contig_13045", "contig_2450", "contig_2464", "contig_2452", "contig_5040",
                         "contig_8600", "contig_7385", "contig_48956"]
    },

    'Caliciviridae': {
        'all_contigs': ["contig_7353", "contig_37552", "contig_6302", "contig_3246", "contig_48995",
                       "contig_48991", "contig_3242", "contig_11605", "contig_3640", "contig_4341",
                       "contig_3213", "contig_22772", "contig_5160", "contig_3210", "contig_3270",
                       "contig_9439", "contig_3231", "contig_3211", "contig_30594", "contig_14281",
                       "contig_7879", "contig_3212", "contig_9179"],
        'US_contigs': ["contig_37552", "contig_3246", "contig_48995", "contig_48991", "contig_3242",
                      "contig_4341", "contig_3270", "contig_9439", "contig_3211", "contig_30594", "contig_3212"],
        'Italy_contigs': ["contig_7353", "contig_6302", "contig_11605", "contig_3640", "contig_4341",
                         "contig_22772", "contig_5160", "contig_3231", "contig_14281", "contig_7879", "contig_9179"]
    },

    'Picornaviridae': {
        'all_contigs': ["contig_3146", "contig_3249", "contig_3176", "contig_292", "contig_49004",
                       "contig_3235", "contig_24559", "contig_10196", "contig_3299", "contig_8104",
                       "contig_33655", "contig_9975", "contig_4215", "contig_3279", "contig_49016",
                       "contig_27255", "contig_43881", "contig_3229", "contig_9150", "contig_3260",
                       "contig_3301", "contig_3290", "contig_13618", "contig_3516", "contig_3096",
                       "contig_3263", "contig_3303", "contig_3225", "contig_44945", "contig_42204",
                       "contig_10454", "contig_49003", "contig_3224", "contig_3287", "contig_3292",
                       "contig_3302", "contig_3280", "contig_18302", "contig_3276", "contig_3241",
                       "contig_3295", "contig_49002", "contig_30058", "contig_3233", "contig_3069",
                       "contig_9184", "contig_11106", "contig_17084", "contig_29259", "contig_3214",
                       "contig_10594", "contig_15445", "contig_37308", "contig_48979", "contig_4110"],
        'US_contigs': ["contig_3146", "contig_3176", "contig_292", "contig_33655", "contig_9975",
                      "contig_4215", "contig_3279", "contig_49016", "contig_27255", "contig_43881",
                      "contig_3229", "contig_9150", "contig_3263", "contig_3303", "contig_3225",
                      "contig_3292", "contig_3302", "contig_18302", "contig_49002", "contig_30058",
                      "contig_3069", "contig_11106", "contig_10594", "contig_48979"],
        'Italy_contigs': ["contig_3176", "contig_49004", "contig_3235", "contig_24559", "contig_10196",
                         "contig_3299", "contig_3260", "contig_13618", "contig_3516", "contig_3096",
                         "contig_44945", "contig_42204", "contig_10454", "contig_49003", "contig_3224",
                         "contig_3287", "contig_3302", "contig_3280", "contig_3276", "contig_3241",
                         "contig_3295", "contig_30058", "contig_3233", "contig_9184", "contig_17084",
                         "contig_29259", "contig_3214", "contig_15445", "contig_37308"]
    }
}

def create_cohort_csv(family_name, family_data, output_dir):
    """
    Create a CSV file for a viral family with cohort information

    Args:
        family_name: Name of the viral family
        family_data: Dictionary with all_contigs, US_contigs, Italy_contigs
        output_dir: Directory to save the CSV file
    """
    all_contigs = family_data['all_contigs']
    us_contigs = set(family_data['US_contigs'])
    italy_contigs = set(family_data['Italy_contigs'])

    # Create dataframe
    data = []
    for contig in all_contigs:
        us_present = "YES" if contig in us_contigs else "NO"
        italy_present = "YES" if contig in italy_contigs else "NO"

        # filtered = YES if both US and Italy are NO, otherwise NO
        filtered = "YES" if (us_present == "NO" and italy_present == "NO") else "NO"

        data.append({
            'Contig': contig,
            'US': us_present,
            'Italy': italy_present,
            'filtered': filtered
        })

    df = pd.DataFrame(data)

    # Save to CSV
    output_file = output_dir / f"{family_name}_contig_cohorts.csv"
    df.to_csv(output_file, index=False)

    print(f"✓ {family_name:15} -> {output_file.name}")
    print(f"  Total contigs: {len(all_contigs)}")
    print(f"  US only: {len([d for d in data if d['US'] == 'YES' and d['Italy'] == 'NO'])}")
    print(f"  Italy only: {len([d for d in data if d['US'] == 'NO' and d['Italy'] == 'YES'])}")
    print(f"  Both cohorts: {len([d for d in data if d['US'] == 'YES' and d['Italy'] == 'YES'])}")
    print(f"  Filtered (neither): {len([d for d in data if d['filtered'] == 'YES'])}")
    print()

def main():
    print("="*70)
    print("Creating contig cohort CSV files for each viral family")
    print("="*70)
    print()

    for family_name, family_data in viral_families.items():
        family_dir = base_dir / family_name
        create_cohort_csv(family_name, family_data, family_dir)

    print("="*70)
    print(f"All {len(viral_families)} CSV files created successfully!")
    print("="*70)

if __name__ == "__main__":
    main()
