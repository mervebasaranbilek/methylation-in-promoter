
import argparse
import csv
import re
from Bio import SeqIO

def is_tss_feature(feat, relaxed=False):
    q = feat.qualifiers

    note_is_tss = any("Transcription Start Site" in s for s in q.get("note", []))
    locus_is_tss = any(str(s).startswith("TSS_") for s in q.get("locus_tag", []))

    base_match = (feat.type == "misc_feature") and (note_is_tss or locus_is_tss)

    if not relaxed:
        return base_match

    type_is_tss = feat.type.lower() in {"tss", "transcription_start_site"}
    return base_match or type_is_tss

def embl_to_csv(in_embl, out_csv, relaxed=False):
    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "tss_position", "strand", "locus_tag", "experiment",
            "citation", "note", "gene", "product", "record_id"
        ])

        for record in SeqIO.parse(in_embl, "embl"):
            for feat in record.features:
                if not is_tss_feature(feat, relaxed=relaxed):
                    continue

                
                pos_1based = int(feat.location.start) + 1
                strand = "+" if feat.location.strand == 1 else "-"

                q = feat.qualifiers
                locus_tag = q.get("locus_tag", [""])[0]
                experiment = q.get("experiment", [""])[0]
                citation = q.get("citation", [""])[0]
                note = q.get("note", [""])[0]
                gene = q.get("gene", [""])[0]
                product = q.get("product", [""])[0]

                writer.writerow([
                    pos_1based, strand, locus_tag, experiment,
                    citation, note, gene, product, record.id
                ])

def main():
    ap = argparse.ArgumentParser(description="Extract TSS entries from EMBL to CSV.")
    ap.add_argument("-i", "--input", required=True, help="Input EMBL file")
    ap.add_argument("-o", "--output", default="tss_data.csv", help="Output CSV")
    ap.add_argument("--relaxed", action="store_true",
                    help="Also include features typed as 'TSS' or 'transcription_start_site'")
    args = ap.parse_args()

    embl_to_csv(args.input, args.output, relaxed=args.relaxed)

if __name__ == "__main__":
    main()
