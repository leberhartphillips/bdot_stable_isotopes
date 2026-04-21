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
    repeat_pair_detected=False,
    repeat_detection_sources="",
    repeat_detection_pattern="",
    repeat_resolution_status="",
    strict_row_action="",
    sensitivity_row_action="",
    sensitivity_selected_source_row="",
    sensitivity_selected_value="",
    audit_retained=True,
    excluded_by_lab_instruction=False,
    contributes_to_model_value=False,
    model_value_method="",
    model_source_rows="",
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
        "repeat_pair_detected": repeat_pair_detected,
        "repeat_detection_sources": repeat_detection_sources,
        "repeat_detection_pattern": repeat_detection_pattern,
        "repeat_resolution_status": repeat_resolution_status,
        "strict_row_action": strict_row_action,
        "sensitivity_row_action": sensitivity_row_action,
        "sensitivity_selected_source_row": sensitivity_selected_source_row,
        "sensitivity_selected_value": (
            sensitivity_selected_value
            if sensitivity_selected_value not in ("", None)
            else ""
        ),
        "audit_retained": audit_retained,
        "excluded_by_lab_instruction": excluded_by_lab_instruction,
        "contributes_to_model_value": contributes_to_model_value,
        "model_value_method": model_value_method,
        "model_source_rows": model_source_rows,
    }


def classify_green_pattern(rows):
    green_flags = [bool(row.get("row_has_green_any")) for row in rows]
    if not any(green_flags):
        return "no_green"
    if all(green_flags):
        return "all_rows_green"

    green_rows = [row for row in rows if row.get("row_has_green_any")]
    if green_rows and all(row.get("is_repeat") for row in green_rows):
        return "rpt_only_green"
    if green_rows and all(not row.get("is_repeat") for row in green_rows):
        return "non_rpt_only_green"
    return "mixed_green"


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
    repeat_pairs = []

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
    repeat_pair_counter = 1

    for sample_base, group_rows in sorted(groups.items()):
        group_rows = sorted(group_rows, key=lambda row: row["source_sheet_row"])
        if len(group_rows) < 2:
            continue

        for row in group_rows:
            row["row_has_green_any"] = any(
                [
                    row["d15n_fill"] == GREEN_FILL_INDEX,
                    row["avg_d15n_fill"] == GREEN_FILL_INDEX and row["avg_d15n"] is not None,
                    row["d13c_fill"] == GREEN_FILL_INDEX,
                    row["avg_d13c_fill"] == GREEN_FILL_INDEX and row["avg_d13c"] is not None,
                ]
            )

        pair_has_rpt_suffix = any(row["is_repeat"] for row in group_rows)
        pair_has_green = any(row["row_has_green_any"] for row in group_rows)
        pair_has_comments = False
        repeat_detection_sources_pair = unique_preserving_order(
            [
                "green formatting" if pair_has_green else "",
                "_rpt suffix" if pair_has_rpt_suffix else "",
            ]
        )
        repeat_detection_pattern = classify_green_pattern(group_rows)
        repeat_detection_pattern_label = (
            "no_green_rpt_named"
            if repeat_detection_pattern == "no_green" and pair_has_rpt_suffix
            else repeat_detection_pattern
        )
        pair_isotopes_present = []
        pair_resolution_statuses = []
        pair_notes = []

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
                explicit_do_not_use = "do not use" in lower_comment
                explicit_use = ("use these" in lower_comment) and (
                    "do not use" not in lower_comment
                )
                pair_has_comments = pair_has_comments or explicit_do_not_use or explicit_use
                annotated_rows.append(
                    {
                        **row,
                        "comment_text": comment_text,
                        "explicit_do_not_use": explicit_do_not_use,
                        "explicit_use": explicit_use,
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
            isotope_detection_by_green = any(
                row["green_raw"] or row["green_avg"] for row in annotated_rows
            )
            isotope_detection_by_comments = any(
                row["explicit_use"] or row["explicit_do_not_use"] for row in annotated_rows
            )
            repeat_pair_detected = (
                len(annotated_rows) >= 2
                and (
                    pair_has_rpt_suffix
                    or isotope_detection_by_green
                    or isotope_detection_by_comments
                )
            )
            repeat_detection_sources = "|".join(
                unique_preserving_order(
                    [
                        "green formatting" if isotope_detection_by_green else "",
                        "_rpt suffix" if pair_has_rpt_suffix else "",
                        "comments" if isotope_detection_by_comments else "",
                    ]
                )
            )

            if not repeat_pair_detected:
                continue

            final_value = None
            final_action = None
            final_reason = None
            evidence_source = None
            rule_basis = None
            notes_text = None
            repeat_resolution_status = None
            strict_model_action = "include_adjudicated"
            sensitivity_model_action = "include_adjudicated"
            sensitivity_selected_source_row = ""
            sensitivity_selected_value = ""

            pair_isotopes_present.append(config["isotope"])

            if len(explicit_keep_rows) > 1 or len(green_average_values) > 1:
                repeat_resolution_status = "conflicting_repeat_instructions"
                final_action = "manual_review"
                final_reason = (
                    "Multiple conflicting workbook instructions were detected for this "
                    "repeat-analysis isotope."
                )
                evidence_source = repeat_detection_sources
                rule_basis = "conflicting_repeat_instructions"
                notes_text = (
                    "Conflicting explicit use instructions or multiple green averages "
                    "require manual review."
                )

                for row in annotated_rows:
                    audit_retained = not row["explicit_do_not_use"]
                    action = "exclude" if row["explicit_do_not_use"] else "keep"
                    reason = (
                        row["comment_text"]
                        or (
                            f"Explicit do-not-use comment for {config['display_name']}."
                            if row["explicit_do_not_use"]
                            else final_reason
                        )
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
                            adjudicated_value="",
                            action=action,
                            reason=reason,
                            evidence_source=evidence_source,
                            decision_class="manual_review",
                            rule_basis=rule_basis,
                            notes=notes_text,
                            repeat_pair_detected=True,
                            repeat_detection_sources=repeat_detection_sources,
                            repeat_detection_pattern=repeat_detection_pattern_label,
                            repeat_resolution_status=repeat_resolution_status,
                            strict_row_action="exclude" if row["explicit_do_not_use"] else "",
                            sensitivity_row_action="exclude" if row["explicit_do_not_use"] else "",
                            sensitivity_selected_source_row="",
                            sensitivity_selected_value="",
                            audit_retained=audit_retained,
                            excluded_by_lab_instruction=row["explicit_do_not_use"],
                            contributes_to_model_value=False,
                            model_value_method="",
                            model_source_rows="",
                        )
                    )
                    decision_counter += 1

                manual_rows.append(
                    {
                        "raw_file": raw_file,
                        "sheet": sheet_name,
                        "sample_identifier": sample_base,
                        "sample_base_identifier": sample_base,
                        "tissue_feather_type": numeric_rows[0]["tissue_feather_type"],
                        "isotope": config["isotope"],
                        "decision_status": "pending_manual_review",
                        "current_build_behavior": (
                            "Retain non-excluded repeat rows in assay-level audit tables, "
                            "but do not derive a modelling value until conflicting workbook "
                            "instructions are resolved."
                        ),
                        "reason": final_reason,
                        "repeat_pair_detected": True,
                        "repeat_detection_sources": repeat_detection_sources,
                        "repeat_detection_pattern": repeat_detection_pattern_label,
                        "strict_model_action": "manual_review_required",
                        "sensitivity_model_action": "manual_review_required",
                        "sensitivity_selected_source_row": "",
                        "sensitivity_selected_value": "",
                        "notes": "Raw source rows: "
                        + "|".join(str(row["source_sheet_row"]) for row in annotated_rows),
                    }
                )

            elif len(explicit_keep_rows) == 1:
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
                    f"{config['isotope']}. Other non-excluded repeat rows are retained "
                    "for assay diagnostics but do not define the modelling value."
                )
                repeat_resolution_status = "resolved_comment_selected_row"
                model_value_method = "explicit_use_selected_row"
                model_source_rows = str(selected_row["source_sheet_row"])

                for row in annotated_rows:
                    if row["explicit_do_not_use"]:
                        action = "exclude"
                        reason = (
                            row["comment_text"]
                            or f"Explicit do-not-use comment for {config['display_name']}."
                        )
                        audit_retained = False
                        contributes_to_model_value = False
                        row_notes = (
                            f"Excluded while source row {selected_row['source_sheet_row']} "
                            "defined the modelling value."
                        )
                    elif row["source_sheet_row"] == selected_row["source_sheet_row"]:
                        action = "keep"
                        reason = final_reason
                        audit_retained = True
                        contributes_to_model_value = True
                        row_notes = notes_text
                    else:
                        action = "keep"
                        reason = (
                            "Retained as a valid repeat observation, but an explicit workbook "
                            "comment selected another row for the modelling value."
                        )
                        audit_retained = True
                        contributes_to_model_value = False
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
                            repeat_pair_detected=True,
                            repeat_detection_sources=repeat_detection_sources,
                            repeat_detection_pattern=repeat_detection_pattern_label,
                            repeat_resolution_status=repeat_resolution_status,
                            strict_row_action="exclude" if action == "exclude" else "include",
                            sensitivity_row_action="exclude" if action == "exclude" else "include",
                            sensitivity_selected_source_row="",
                            sensitivity_selected_value="",
                            audit_retained=audit_retained,
                            excluded_by_lab_instruction=row["explicit_do_not_use"],
                            contributes_to_model_value=contributes_to_model_value,
                            model_value_method=model_value_method,
                            model_source_rows=model_source_rows,
                        )
                    )
                    decision_counter += 1

            elif len(green_average_values) == 1 and len(usable_rows) >= 1:
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
                    f"Sheet note: {sheet_note_green} "
                    "Non-excluded repeat rows are retained for assay diagnostics."
                )
                repeat_resolution_status = "resolved_use_average"
                model_value_method = "lab_green_average"
                model_source_rows = "|".join(
                    str(row["source_sheet_row"]) for row in usable_rows
                )

                for row in annotated_rows:
                    if row["explicit_do_not_use"]:
                        action = "exclude"
                        reason = (
                            row["comment_text"]
                            or f"Explicit do-not-use comment for {config['display_name']}."
                        )
                        audit_retained = False
                        contributes_to_model_value = False
                        row_notes = (
                            f"Excluded before applying green average from row "
                            f"{average_source_row}."
                        )
                    else:
                        action = "keep"
                        reason = final_reason
                        audit_retained = True
                        contributes_to_model_value = True
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
                            repeat_pair_detected=True,
                            repeat_detection_sources=repeat_detection_sources,
                            repeat_detection_pattern=repeat_detection_pattern_label,
                            repeat_resolution_status=repeat_resolution_status,
                            strict_row_action="exclude" if action == "exclude" else "include",
                            sensitivity_row_action="exclude" if action == "exclude" else "include",
                            sensitivity_selected_source_row="",
                            sensitivity_selected_value="",
                            audit_retained=audit_retained,
                            excluded_by_lab_instruction=row["explicit_do_not_use"],
                            contributes_to_model_value=contributes_to_model_value,
                            model_value_method=model_value_method,
                            model_source_rows=model_source_rows,
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
                repeat_resolution_status = "resolved_comment_single_remaining"
                model_value_method = "single_remaining_after_exclusion"
                model_source_rows = str(selected_row["source_sheet_row"])

                for row in annotated_rows:
                    if row["source_sheet_row"] == selected_row["source_sheet_row"]:
                        action = "keep"
                        reason = final_reason
                        audit_retained = True
                        contributes_to_model_value = True
                        row_notes = notes_text
                    else:
                        action = "exclude"
                        reason = (
                            row["comment_text"]
                            or f"Explicit do-not-use comment for {config['display_name']}."
                        )
                        audit_retained = False
                        contributes_to_model_value = False
                        row_notes = (
                            f"Excluded while source row {selected_row['source_sheet_row']} "
                            "supplied the modelling value."
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
                            repeat_pair_detected=True,
                            repeat_detection_sources=repeat_detection_sources,
                            repeat_detection_pattern=repeat_detection_pattern_label,
                            repeat_resolution_status=repeat_resolution_status,
                            strict_row_action="exclude" if action == "exclude" else "include",
                            sensitivity_row_action="exclude" if action == "exclude" else "include",
                            sensitivity_selected_source_row="",
                            sensitivity_selected_value="",
                            audit_retained=audit_retained,
                            excluded_by_lab_instruction=row["explicit_do_not_use"],
                            contributes_to_model_value=contributes_to_model_value,
                            model_value_method=model_value_method,
                            model_source_rows=model_source_rows,
                        )
                    )
                    decision_counter += 1

            elif len(usable_rows) >= 1:
                final_value = sum(row[config["raw_key"]] for row in usable_rows) / len(
                    usable_rows
                )
                final_action = "keep"
                final_reason = (
                    "No explicit workbook instruction selected one repeat analysis, so "
                    "the modelling value is the mean of all retained valid repeat assays."
                )
                evidence_source = repeat_detection_sources
                rule_basis = "mean_of_valid_repeats"
                notes_text = (
                    "Retained all non-excluded repeat rows for assay diagnostics and "
                    "used their mean as the adjudicated modelling value."
                )
                repeat_resolution_status = "resolved_mean_of_valid_repeats"
                model_value_method = "mean_of_valid_repeats"
                model_source_rows = "|".join(
                    str(row["source_sheet_row"]) for row in usable_rows
                )

                for row in annotated_rows:
                    if row["explicit_do_not_use"]:
                        action = "exclude"
                        reason = (
                            row["comment_text"]
                            or f"Explicit do-not-use comment for {config['display_name']}."
                        )
                        audit_retained = False
                        contributes_to_model_value = False
                        row_notes = (
                            "Excluded by explicit workbook instruction before calculating "
                            "the mean of valid repeats."
                        )
                    else:
                        action = "keep"
                        reason = final_reason
                        audit_retained = True
                        contributes_to_model_value = True
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
                            repeat_pair_detected=True,
                            repeat_detection_sources=repeat_detection_sources,
                            repeat_detection_pattern=repeat_detection_pattern_label,
                            repeat_resolution_status=repeat_resolution_status,
                            strict_row_action="exclude" if action == "exclude" else "include",
                            sensitivity_row_action="exclude" if action == "exclude" else "include",
                            sensitivity_selected_source_row="",
                            sensitivity_selected_value="",
                            audit_retained=audit_retained,
                            excluded_by_lab_instruction=row["explicit_do_not_use"],
                            contributes_to_model_value=contributes_to_model_value,
                            model_value_method=model_value_method,
                            model_source_rows=model_source_rows,
                        )
                    )
                    decision_counter += 1

            else:
                repeat_resolution_status = "all_rows_explicitly_excluded"
                final_action = "exclude"
                final_reason = (
                    "All repeat rows for this isotope were explicitly excluded by "
                    "workbook instructions."
                )
                evidence_source = "comment column"
                rule_basis = "all_repeat_rows_explicitly_excluded"
                notes_text = (
                    "No retained repeat rows remain for assay diagnostics or modelling."
                )

                for row in annotated_rows:
                    action = "exclude"
                    reason = (
                        row["comment_text"]
                        or f"Explicit do-not-use comment for {config['display_name']}."
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
                            adjudicated_value="",
                            action=action,
                            reason=reason,
                            evidence_source=evidence_source,
                            decision_class="rule_based",
                            rule_basis=rule_basis,
                            notes=notes_text,
                            repeat_pair_detected=True,
                            repeat_detection_sources=repeat_detection_sources,
                            repeat_detection_pattern=repeat_detection_pattern_label,
                            repeat_resolution_status=repeat_resolution_status,
                            strict_row_action="exclude",
                            sensitivity_row_action="exclude",
                            sensitivity_selected_source_row="",
                            sensitivity_selected_value="",
                            audit_retained=False,
                            excluded_by_lab_instruction=True,
                            contributes_to_model_value=False,
                            model_value_method="",
                            model_source_rows="",
                        )
                    )
                    decision_counter += 1

            pair_resolution_statuses.append(repeat_resolution_status)
            pair_notes.append(notes_text)

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
                    "repeat_pair_detected": True,
                    "repeat_detection_sources": repeat_detection_sources,
                    "repeat_detection_pattern": repeat_detection_pattern_label,
                    "repeat_resolution_status": repeat_resolution_status,
                    "strict_model_action": strict_model_action,
                    "sensitivity_model_action": sensitivity_model_action,
                    "sensitivity_selected_source_row": sensitivity_selected_source_row,
                    "sensitivity_selected_value": sensitivity_selected_value,
                    "model_value_method": rule_basis if final_value is not None else "",
                }
            )

        overall_sources = unique_preserving_order(
            repeat_detection_sources_pair
            + (["comments"] if pair_has_comments else [])
        )
        repeat_pairs.append(
            {
                "repeat_pair_id": f"repeat_pair_{repeat_pair_counter:03d}",
                "raw_file": raw_file,
                "sheet": sheet_name,
                "sample_identifier": sample_base,
                "sample_base_identifier": sample_base,
                "tissue_feather_type": group_rows[0]["tissue_feather_type"],
                "n_rows_in_pair": len(group_rows),
                "raw_source_sheet_rows": "|".join(
                    str(row["source_sheet_row"]) for row in group_rows
                ),
                "sample_identifiers_raw": "|".join(
                    row["sample_identifier"] for row in group_rows
                ),
                "repeat_pair_detected": True,
                "detection_by_green_formatting": pair_has_green,
                "detection_by_rpt_suffix": pair_has_rpt_suffix,
                "detection_by_comments": pair_has_comments,
                "repeat_detection_sources": "|".join(overall_sources),
                "repeat_detection_pattern": repeat_detection_pattern_label,
                "isotopes_present": "|".join(unique_preserving_order(pair_isotopes_present)),
                "resolution_statuses": "|".join(unique_preserving_order(pair_resolution_statuses)),
                "notes": " | ".join(unique_preserving_order(pair_notes)),
            }
        )
        repeat_pair_counter += 1

    repeat_pair_summary = []
    if repeat_pairs:
        repeat_pair_summary.extend(
            [
                {
                    "metric": "detected_by_green_formatting",
                    "value": "true",
                    "n_repeat_pairs": sum(
                        1
                        for row in repeat_pairs
                        if row["detection_by_green_formatting"]
                    ),
                },
                {
                    "metric": "detected_by_rpt_suffix",
                    "value": "true",
                    "n_repeat_pairs": sum(
                        1 for row in repeat_pairs if row["detection_by_rpt_suffix"]
                    ),
                },
                {
                    "metric": "detected_by_comments",
                    "value": "true",
                    "n_repeat_pairs": sum(
                        1 for row in repeat_pairs if row["detection_by_comments"]
                    ),
                },
            ]
        )
        pattern_counts = defaultdict(int)
        for row in repeat_pairs:
            pattern_counts[row["repeat_detection_pattern"]] += 1
        for pattern, count in sorted(pattern_counts.items()):
            repeat_pair_summary.append(
                {
                    "metric": "repeat_detection_pattern",
                    "value": pattern,
                    "n_repeat_pairs": count,
                }
            )

    return decision_rows, affected_samples, manual_rows, repeat_pairs, repeat_pair_summary


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
                    audit_retained=False,
                    excluded_by_lab_instruction=True,
                    contributes_to_model_value=False,
                    model_value_method="",
                    model_source_rows="",
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
    (
        live_decisions,
        live_affected,
        manual_rows,
        repeat_pairs,
        repeat_pair_summary,
    ) = parse_live_batch5_cn()
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
        "repeat_pair_detected",
        "repeat_detection_sources",
        "repeat_detection_pattern",
        "repeat_resolution_status",
        "strict_row_action",
        "sensitivity_row_action",
        "sensitivity_selected_source_row",
        "sensitivity_selected_value",
        "audit_retained",
        "excluded_by_lab_instruction",
        "contributes_to_model_value",
        "model_value_method",
        "model_source_rows",
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
        "repeat_pair_detected",
        "repeat_detection_sources",
        "repeat_detection_pattern",
        "repeat_resolution_status",
        "strict_model_action",
        "sensitivity_model_action",
        "sensitivity_selected_source_row",
        "sensitivity_selected_value",
        "model_value_method",
    ]
    manual_fieldnames = [
        "raw_file",
        "sheet",
        "sample_identifier",
        "sample_base_identifier",
        "tissue_feather_type",
        "isotope",
        "decision_status",
        "current_build_behavior",
        "reason",
        "repeat_pair_detected",
        "repeat_detection_sources",
        "repeat_detection_pattern",
        "strict_model_action",
        "sensitivity_model_action",
        "sensitivity_selected_source_row",
        "sensitivity_selected_value",
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
    repeat_pair_fieldnames = [
        "repeat_pair_id",
        "raw_file",
        "sheet",
        "sample_identifier",
        "sample_base_identifier",
        "tissue_feather_type",
        "n_rows_in_pair",
        "raw_source_sheet_rows",
        "sample_identifiers_raw",
        "repeat_pair_detected",
        "detection_by_green_formatting",
        "detection_by_rpt_suffix",
        "detection_by_comments",
        "repeat_detection_sources",
        "repeat_detection_pattern",
        "isotopes_present",
        "resolution_statuses",
        "notes",
    ]
    repeat_pair_summary_fieldnames = [
        "metric",
        "value",
        "n_repeat_pairs",
    ]

    write_csv(DERIVED_DIR / "assay_qc_adjudication.csv", adjudication_rows, adjudication_fieldnames)
    write_csv(DERIVED_DIR / "assay_qc_affected_samples.csv", affected_rows, affected_fieldnames)
    write_csv(DERIVED_DIR / "assay_qc_manual_adjudication.csv", manual_rows, manual_fieldnames)
    write_csv(DERIVED_DIR / "assay_qc_sheet_rules.csv", build_sheet_rules(), rule_fieldnames)
    write_csv(DERIVED_DIR / "assay_qc_repeat_pairs.csv", repeat_pairs, repeat_pair_fieldnames)
    write_csv(
        DERIVED_DIR / "assay_qc_repeat_pair_summary.csv",
        repeat_pair_summary,
        repeat_pair_summary_fieldnames,
    )


if __name__ == "__main__":
    main()
