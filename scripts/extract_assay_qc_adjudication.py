#!/usr/bin/env python3

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "_vendor_py"))

import xlrd  # noqa: E402

DERIVED_DIR = REPO_ROOT / "data" / "derived"

GREEN_FILL_INDEX = 42
RED_FONT_INDEX = 10


def normalize_space(value):
    if value is None:
        return ""
    text = str(value).replace("\xa0", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def normalize_lower(value):
    return normalize_space(value).lower()


def numeric_or_none(value):
    if value in ("", None):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def unique_preserving_order(values):
    seen = set()
    out = []
    for value in values:
        if value in ("", None):
            continue
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def standardize_tissue(label):
    label = normalize_lower(label)
    if "breast" in label:
        return "Breast"
    if "primary" in label or "pimary" in label:
        return "Primary"
    return ""


def strip_repeat_suffix(sample_identifier):
    sample_identifier = normalize_space(sample_identifier)
    return re.sub(r"_rpt$", "", sample_identifier, flags=re.IGNORECASE)


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def get_style(workbook, sheet, rowx, colx):
    xf = workbook.xf_list[sheet.cell_xf_index(rowx, colx)]
    font = workbook.font_list[xf.font_index]
    return {
        "fill_index": xf.background.pattern_colour_index,
        "font_index": font.colour_index,
    }


def build_decision_row(
    *,
    decision_id,
    raw_file,
    sheet,
    source_sheet_row,
    sample_identifier,
    sample_base_identifier,
    tissue_feather_type,
    isotope,
    raw_value,
    adjudicated_value,
    action,
    reason,
    evidence_source,
    decision_class,
    rule_basis,
    notes,
):
    return {
        "decision_id": decision_id,
        "raw_file": raw_file,
        "sheet": sheet,
        "source_sheet_row": source_sheet_row,
        "sample_identifier": sample_identifier,
        "sample_base_identifier": sample_base_identifier,
        "tissue_feather_type": tissue_feather_type,
        "isotope": isotope,
        "raw_value": raw_value if raw_value is not None else "",
        "adjudicated_value": adjudicated_value if adjudicated_value is not None else "",
        "action": action,
        "reason": reason,
        "evidence_source": evidence_source,
        "decision_class": decision_class,
        "rule_basis": rule_basis,
        "notes": notes,
    }


def parse_live_batch5_cn():
    raw_file = "data/Master_EA_C_N_SI_Batch 5 Breast and Primary Feathers_LEH_Repeats added_April25_v1_SB_v2.xls"
    sheet_name = "Sample results_C&N"
    workbook = xlrd.open_workbook(REPO_ROOT / raw_file, formatting_info=True)
    sheet = workbook.sheet_by_name(sheet_name)

    notes = [
        normalize_space(sheet.cell_value(0, 1)),
        normalize_space(sheet.cell_value(5, 6)),
    ]
    sheet_note_green = notes[0]

    rows = []
    for rowx in range(12, sheet.nrows):
        sample_identifier = normalize_space(sheet.cell_value(rowx, 1))
        if not sample_identifier:
            continue

        d15n_style = get_style(workbook, sheet, rowx, 12)
        avg_d15n_style = get_style(workbook, sheet, rowx, 15)
        d13c_style = get_style(workbook, sheet, rowx, 22)
        avg_d13c_style = get_style(workbook, sheet, rowx, 26)

        rows.append(
            {
                "source_sheet_row": rowx,
                "sample_identifier": sample_identifier,
                "sample_base_identifier": strip_repeat_suffix(sample_identifier),
                "tissue_feather_type": standardize_tissue(sample_identifier),
                "is_repeat": bool(normalize_space(sheet.cell_value(rowx, 2)))
                or sample_identifier.lower().endswith("_rpt"),
                "d15n": numeric_or_none(sheet.cell_value(rowx, 12)),
                "avg_d15n": numeric_or_none(sheet.cell_value(rowx, 15)),
                "d15n_fill": d15n_style["fill_index"],
                "d15n_font": d15n_style["font_index"],
                "avg_d15n_fill": avg_d15n_style["fill_index"],
                "n_comment": normalize_space(sheet.cell_value(rowx, 17)),
                "d13c": numeric_or_none(sheet.cell_value(rowx, 22)),
                "avg_d13c": numeric_or_none(sheet.cell_value(rowx, 26)),
                "d13c_fill": d13c_style["fill_index"],
                "d13c_font": d13c_style["font_index"],
                "avg_d13c_fill": avg_d13c_style["fill_index"],
                "c_comment": normalize_space(sheet.cell_value(rowx, 30)),
                "c_comment_aux": normalize_space(sheet.cell_value(rowx, 28)),
            }
        )

    groups = defaultdict(list)
    for row in rows:
        groups[row["sample_base_identifier"]].append(row)

    decision_rows = []
    affected_samples = []
    manual_rows = []

    isotope_configs = [
        {
            "isotope": "normalised_d15n",
            "raw_key": "d15n",
            "avg_key": "avg_d15n",
            "fill_key": "d15n_fill",
            "font_key": "d15n_font",
            "comment_keys": ["n_comment"],
            "display_name": "N",
        },
        {
            "isotope": "normalised_d13c",
            "raw_key": "d13c",
            "avg_key": "avg_d13c",
            "fill_key": "d13c_fill",
            "font_key": "d13c_font",
            "comment_keys": ["c_comment", "c_comment_aux"],
            "display_name": "C",
        },
    ]

    decision_counter = 1

    for sample_base, group_rows in sorted(groups.items()):
        group_rows = sorted(group_rows, key=lambda row: row["source_sheet_row"])
        if len(group_rows) < 2 or not any(row["is_repeat"] for row in group_rows):
            continue

        for config in isotope_configs:
            numeric_rows = [
                row for row in group_rows if row[config["raw_key"]] is not None
            ]
            if len(numeric_rows) < 2:
                continue

            annotated_rows = []
            for row in numeric_rows:
                comment_texts = unique_preserving_order(
                    [row[key] for key in config["comment_keys"]]
                )
                comment_text = " | ".join(comment_texts)
                lower_comment = comment_text.lower()
                annotated_rows.append(
                    {
                        **row,
                        "comment_text": comment_text,
                        "explicit_do_not_use": "do not use" in lower_comment,
                        "explicit_use": ("use these" in lower_comment)
                        and ("do not use" not in lower_comment),
                        "green_raw": row[config["fill_key"]] == GREEN_FILL_INDEX,
                        "red_font_raw": row[config["font_key"]] == RED_FONT_INDEX,
                        "green_avg": row[config["avg_key"]] is not None
                        and row[f"{config['avg_key']}_fill"] == GREEN_FILL_INDEX,
                    }
                )

            explicit_keep_rows = [
                row for row in annotated_rows if row["explicit_use"]
            ]
            explicit_exclude_rows = [
                row for row in annotated_rows if row["explicit_do_not_use"]
            ]
            usable_rows = [
                row for row in annotated_rows if not row["explicit_do_not_use"]
            ]
            green_average_rows = [
                row for row in annotated_rows if row["green_avg"]
            ]
            green_average_values = unique_preserving_order(
                [row[config["avg_key"]] for row in green_average_rows]
            )

            final_value = None
            final_action = None
            final_reason = None
            evidence_source = None
            rule_basis = None
            notes_text = None

            if len(explicit_keep_rows) == 1:
                selected_row = explicit_keep_rows[0]
                final_value = selected_row[config["raw_key"]]
                final_action = "keep"
                final_reason = (
                    f"Repeat-analysis comment explicitly instructs use of the "
                    f"{config['display_name']} value from this assay row."
                )
                evidence_parts = ["comment column"]
                if selected_row["green_raw"]:
                    evidence_parts.append("green raw")
                evidence_source = "|".join(evidence_parts)
                rule_basis = "comment_use_these"
                notes_text = (
                    f"Selected source row {selected_row['source_sheet_row']} for "
                    f"{config['isotope']}."
                )

                for row in annotated_rows:
                    if row["source_sheet_row"] == selected_row["source_sheet_row"]:
                        action = "keep"
                        reason = final_reason
                        row_notes = notes_text
                    elif row["explicit_do_not_use"]:
                        action = "exclude"
                        reason = (
                            row["comment_text"]
                            or f"Explicit do-not-use comment for {config['display_name']}."
                        )
                        row_notes = (
                            f"Excluded while {selected_row['source_sheet_row']} supplied the "
                            "adjudicated value."
                        )
                    else:
                        action = "exclude"
                        reason = (
                            "Superseded by the repeat-analysis row explicitly marked "
                            "for use."
                        )
                        row_notes = (
                            f"Superseded by source row {selected_row['source_sheet_row']}."
                        )

                    decision_rows.append(
                        build_decision_row(
                            decision_id=f"aqc_{decision_counter:04d}",
                            raw_file=raw_file,
                            sheet=sheet_name,
                            source_sheet_row=row["source_sheet_row"],
                            sample_identifier=row["sample_identifier"],
                            sample_base_identifier=sample_base,
                            tissue_feather_type=row["tissue_feather_type"],
                            isotope=config["isotope"],
                            raw_value=row[config["raw_key"]],
                            adjudicated_value=final_value,
                            action=action,
                            reason=reason,
                            evidence_source=evidence_source,
                            decision_class="rule_based",
                            rule_basis=rule_basis,
                            notes=row_notes,
                        )
                    )
                    decision_counter += 1

            elif len(green_average_values) == 1 and len(usable_rows) >= 2:
                final_value = green_average_values[0]
                average_source_row = green_average_rows[0]["source_sheet_row"]
                final_action = "replace_with_average"
                final_reason = (
                    "Green-highlighted average is the lab-approved adjudicated "
                    f"{config['display_name']} value for this repeated assay."
                )
                evidence_source = "sheet note|green average"
                rule_basis = "sheet_note_green_average"
                notes_text = (
                    f"Average stored on source row {average_source_row}. "
                    f"Sheet note: {sheet_note_green}"
                )

                for row in annotated_rows:
                    if row["explicit_do_not_use"]:
                        action = "exclude"
                        reason = (
                            row["comment_text"]
                            or f"Explicit do-not-use comment for {config['display_name']}."
                        )
                        row_notes = (
                            f"Excluded before applying green average from row "
                            f"{average_source_row}."
                        )
                    else:
                        action = "replace_with_average"
                        reason = final_reason
                        row_notes = notes_text

                    decision_rows.append(
                        build_decision_row(
                            decision_id=f"aqc_{decision_counter:04d}",
                            raw_file=raw_file,
                            sheet=sheet_name,
                            source_sheet_row=row["source_sheet_row"],
                            sample_identifier=row["sample_identifier"],
                            sample_base_identifier=sample_base,
                            tissue_feather_type=row["tissue_feather_type"],
                            isotope=config["isotope"],
                            raw_value=row[config["raw_key"]],
                            adjudicated_value=final_value,
                            action=action,
                            reason=reason,
                            evidence_source=evidence_source,
                            decision_class="rule_based",
                            rule_basis=rule_basis,
                            notes=row_notes,
                        )
                    )
                    decision_counter += 1

            elif len(explicit_exclude_rows) >= 1 and len(usable_rows) == 1:
                selected_row = usable_rows[0]
                final_value = selected_row[config["raw_key"]]
                final_action = "keep"
                final_reason = (
                    "Only one usable assay remained after explicit do-not-use "
                    "exclusions."
                )
                evidence_source = "comment column"
                rule_basis = "comment_do_not_use_single_remaining"
                notes_text = (
                    f"Selected remaining usable source row {selected_row['source_sheet_row']}."
                )

                for row in annotated_rows:
                    if row["source_sheet_row"] == selected_row["source_sheet_row"]:
                        action = "keep"
                        reason = final_reason
                        row_notes = notes_text
                    else:
                        action = "exclude"
                        reason = (
                            row["comment_text"]
                            or f"Explicit do-not-use comment for {config['display_name']}."
                        )
                        row_notes = (
                            f"Excluded while {selected_row['source_sheet_row']} supplied the "
                            "adjudicated value."
                        )

                    decision_rows.append(
                        build_decision_row(
                            decision_id=f"aqc_{decision_counter:04d}",
                            raw_file=raw_file,
                            sheet=sheet_name,
                            source_sheet_row=row["source_sheet_row"],
                            sample_identifier=row["sample_identifier"],
                            sample_base_identifier=sample_base,
                            tissue_feather_type=row["tissue_feather_type"],
                            isotope=config["isotope"],
                            raw_value=row[config["raw_key"]],
                            adjudicated_value=final_value,
                            action=action,
                            reason=reason,
                            evidence_source=evidence_source,
                            decision_class="rule_based",
                            rule_basis=rule_basis,
                            notes=row_notes,
                        )
                    )
                    decision_counter += 1

            else:
                manual_rows.append(
                    {
                        "raw_file": raw_file,
                        "sheet": sheet_name,
                        "sample_identifier": sample_base,
                        "tissue_feather_type": numeric_rows[0]["tissue_feather_type"],
                        "isotope": config["isotope"],
                        "decision_status": "pending_manual_review",
                        "current_build_behavior": (
                            "No assay-level override applied; retain raw usable assay "
                            "values separately and flag downstream aggregation as unresolved."
                        ),
                        "reason": (
                            "Repeated assay present but no parseable green average or "
                            "explicit use/exclude instruction was found for this isotope."
                        ),
                        "notes": (
                            "Raw source rows: "
                            + "|".join(str(row["source_sheet_row"]) for row in annotated_rows)
                        ),
                    }
                )
                continue

            affected_samples.append(
                {
                    "raw_file": raw_file,
                    "sheet": sheet_name,
                    "sample_identifier": sample_base,
                    "sample_base_identifier": sample_base,
                    "tissue_feather_type": numeric_rows[0]["tissue_feather_type"],
                    "isotope": config["isotope"],
                    "raw_source_sheet_rows": "|".join(
                        str(row["source_sheet_row"]) for row in annotated_rows
                    ),
                    "raw_values": "|".join(
                        str(row[config["raw_key"]]) for row in annotated_rows
                    ),
                    "adjudicated_value": final_value if final_value is not None else "",
                    "final_action": final_action,
                    "reason": final_reason,
                    "evidence_source": evidence_source,
                    "decision_class": "rule_based",
                    "rule_basis": rule_basis,
                    "notes": notes_text,
                }
            )

    return decision_rows, affected_samples, manual_rows


def parse_red_font_exclusions():
    decision_rows = []
    affected_samples = []
    decision_counter = 1

    configs = [
        {
            "raw_file": "data/Master_TCEA_O_SI_Breast Feathers_Replication_Batch 4_LEH_Sept25_SB.xls",
            "sheet": "Sample results_O",
            "header_row": 6,
            "data_start_row": 7,
            "sample_col": 1,
            "tissue": "Breast",
            "isotope": "normalised_d18o",
            "value_col": 12,
            "comment_col": 16,
            "sheet_note_row": 4,
            "sheet_note_col": 1,
        },
        {
            "raw_file": "data/Master_TCEA_H_SI_Batch 1_triplicate breast feathers_LEH_Nov25_SB.xls",
            "sheet": "Sample results_H",
            "header_row": 9,
            "data_start_row": 10,
            "sample_col": 1,
            "tissue": "Breast",
            "isotope": "normalised_d2h",
            "value_col": 10,
            "comment_col": 13,
            "sheet_note_row": 5,
            "sheet_note_col": 1,
        },
    ]

    for config in configs:
        workbook = xlrd.open_workbook(REPO_ROOT / config["raw_file"], formatting_info=True)
        sheet = workbook.sheet_by_name(config["sheet"])
        sheet_note = normalize_space(
            sheet.cell_value(config["sheet_note_row"], config["sheet_note_col"])
        )

        for rowx in range(config["data_start_row"], sheet.nrows):
            sample_identifier = normalize_space(sheet.cell_value(rowx, config["sample_col"]))
            raw_value = numeric_or_none(sheet.cell_value(rowx, config["value_col"]))
            if not sample_identifier or raw_value is None:
                continue

            comment_text = normalize_space(sheet.cell_value(rowx, config["comment_col"]))
            style = get_style(workbook, sheet, rowx, config["value_col"])
            has_red_font = style["font_index"] == RED_FONT_INDEX
            has_do_not_use = "do not use" in comment_text.lower()

            if not has_red_font and not has_do_not_use:
                continue

            evidence_parts = []
            if sheet_note:
                evidence_parts.append("sheet note")
            if has_red_font:
                evidence_parts.append("red font")
            if comment_text:
                evidence_parts.append("comment column")

            reason = comment_text or "Red-font isotope value excluded per sheet note."
            notes = sheet_note if sheet_note else ""

            decision_rows.append(
                build_decision_row(
                    decision_id=f"aqc_rule_{decision_counter:04d}",
                    raw_file=config["raw_file"],
                    sheet=config["sheet"],
                    source_sheet_row=rowx + 1,
                    sample_identifier=sample_identifier,
                    sample_base_identifier=sample_identifier,
                    tissue_feather_type=config["tissue"],
                    isotope=config["isotope"],
                    raw_value=raw_value,
                    adjudicated_value=None,
                    action="exclude",
                    reason=reason,
                    evidence_source="|".join(evidence_parts),
                    decision_class="rule_based",
                    rule_basis="sheet_note_red_font_exclusion",
                    notes=notes,
                )
            )
            decision_counter += 1

            affected_samples.append(
                {
                    "raw_file": config["raw_file"],
                    "sheet": config["sheet"],
                    "sample_identifier": sample_identifier,
                    "sample_base_identifier": sample_identifier,
                    "tissue_feather_type": config["tissue"],
                    "isotope": config["isotope"],
                    "raw_source_sheet_rows": str(rowx + 1),
                    "raw_values": str(raw_value),
                    "adjudicated_value": "",
                    "final_action": "exclude",
                    "reason": reason,
                    "evidence_source": "|".join(evidence_parts),
                    "decision_class": "rule_based",
                    "rule_basis": "sheet_note_red_font_exclusion",
                    "notes": notes,
                }
            )

    return decision_rows, affected_samples


def build_sheet_rules():
    return [
        {
            "raw_file": "data/Master_EA_C_N_SI_Batch 5 Breast and Primary Feathers_LEH_Repeats added_April25_v1_SB_v2.xls",
            "sheet": "Sample results_C&N",
            "note_scope": "sheet note",
            "note_text": (
                "10 samples were repeat analysed due to autosampler jamming. "
                "You can use the C and N data of the individual analyses and/or "
                "averages, where data are highlighted in green."
            ),
            "parse_status": "parsed_rule_based",
            "applied_to_build": "yes",
            "notes": "Used to detect lab-approved green averages in repeated batch 5 assays.",
        },
        {
            "raw_file": "data/Master_EA_C_N_SI_Batch 5 Breast and Primary Feathers_LEH_Repeats added_April25_v1_SB_v2.xls",
            "sheet": "Sample results_C&N",
            "note_scope": "sheet note",
            "note_text": "Need to repeat the analysis of one feather listed in green boxes at no extra charge.",
            "parse_status": "parsed_metadata_only",
            "applied_to_build": "no",
            "notes": "Recorded as inventory context only; no adjudicated value exists yet for the missing repeat.",
        },
        {
            "raw_file": "data/Master_TCEA_O_SI_Breast Feathers_Replication_Batch 4_LEH_Sept25_SB.xls",
            "sheet": "Sample results_O",
            "note_scope": "sheet note",
            "note_text": "do not use values in red font",
            "parse_status": "parsed_rule_based",
            "applied_to_build": "yes",
            "notes": "Used to exclude red-font oxygen values, including AV21986_A feather.",
        },
        {
            "raw_file": "data/Master_TCEA_H_SI_Batch 1_triplicate breast feathers_LEH_Nov25_SB.xls",
            "sheet": "Sample results_H",
            "note_scope": "sheet note",
            "note_text": "do not use values in red font",
            "parse_status": "parsed_rule_based",
            "applied_to_build": "yes",
            "notes": "Used to exclude red-font hydrogen values, including AV705_D feather.",
        },
        {
            "raw_file": "data/Master_EA_C_N_SI_Batch 3 vane_barb and Batch 4 Breast Feathers_LEH_Sept25_v1_SB.xls",
            "sheet": "Sample results_C&N",
            "note_scope": "sheet note",
            "note_text": "Samples in red text were on the inventory but not analysed for C&N.",
            "parse_status": "parsed_metadata_only",
            "applied_to_build": "no",
            "notes": "Not used as a value-level override because these rows already lack assay values.",
        },
        {
            "raw_file": "data/Master_EA_C_N_SI_Batch 3 vane_barb and Batch 4 Breast Feathers_LEH_Sept25_v1_SB.xls",
            "sheet": "Sample results_C&N",
            "note_scope": "sheet note",
            "note_text": "Samples in yellow highlight were analysed for shaft-barb data comparison.",
            "parse_status": "parsed_metadata_only",
            "applied_to_build": "no",
            "notes": "Useful provenance metadata, but not a usable-value exclusion rule.",
        },
        {
            "raw_file": "data/old file versions/Master_SI_Feathers_LEH_Sept24_v1_SB_v3_sample wts added.xls",
            "sheet": "Sample results_C&N",
            "note_scope": "sheet note",
            "note_text": (
                "Triplicate note suggests homogenising all 3 feathers together for "
                "future analyses, while preserving context about between-feather variability."
            ),
            "parse_status": "context_only",
            "applied_to_build": "no",
            "notes": "Contextual interpretation only; no explicit per-value keep/exclude instruction was parsed.",
        },
    ]


def main():
    live_decisions, live_affected, manual_rows = parse_live_batch5_cn()
    red_font_decisions, red_font_affected = parse_red_font_exclusions()

    adjudication_rows = live_decisions + red_font_decisions
    affected_rows = live_affected + red_font_affected

    adjudication_fieldnames = [
        "decision_id",
        "raw_file",
        "sheet",
        "source_sheet_row",
        "sample_identifier",
        "sample_base_identifier",
        "tissue_feather_type",
        "isotope",
        "raw_value",
        "adjudicated_value",
        "action",
        "reason",
        "evidence_source",
        "decision_class",
        "rule_basis",
        "notes",
    ]
    affected_fieldnames = [
        "raw_file",
        "sheet",
        "sample_identifier",
        "sample_base_identifier",
        "tissue_feather_type",
        "isotope",
        "raw_source_sheet_rows",
        "raw_values",
        "adjudicated_value",
        "final_action",
        "reason",
        "evidence_source",
        "decision_class",
        "rule_basis",
        "notes",
    ]
    manual_fieldnames = [
        "raw_file",
        "sheet",
        "sample_identifier",
        "tissue_feather_type",
        "isotope",
        "decision_status",
        "current_build_behavior",
        "reason",
        "notes",
    ]
    rule_fieldnames = [
        "raw_file",
        "sheet",
        "note_scope",
        "note_text",
        "parse_status",
        "applied_to_build",
        "notes",
    ]

    write_csv(DERIVED_DIR / "assay_qc_adjudication.csv", adjudication_rows, adjudication_fieldnames)
    write_csv(DERIVED_DIR / "assay_qc_affected_samples.csv", affected_rows, affected_fieldnames)
    write_csv(DERIVED_DIR / "assay_qc_manual_adjudication.csv", manual_rows, manual_fieldnames)
    write_csv(DERIVED_DIR / "assay_qc_sheet_rules.csv", build_sheet_rules(), rule_fieldnames)


if __name__ == "__main__":
    main()
