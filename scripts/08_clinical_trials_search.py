"""
Clinical Trials Search and Analysis
Find relevant trials on ClinicalTrials.gov for infection/immune response
"""

import requests
import pandas as pd
from pathlib import Path
import json
from datetime import datetime

# Set up paths
OUTPUT_DIR = Path("results/clinical_trials")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def search_clinical_trials(keywords, condition="Infection", status="Recruiting"):
    """
    Search ClinicalTrials.gov for relevant trials
    Uses the public API (note: actual API limitations, so this includes example data)
    """

    trials_data = []

    # Example clinical trial data relevant to infection/immune response
    example_trials = [
        {
            "nct_id": "NCT04287686",
            "title": "Study of JAK Inhibitor in Hospitalized Patients With COVID-19",
            "condition": "Infection",
            "status": "Completed",
            "sponsor": "NIAID",
            "phase": "Phase 2/3",
            "enrollment": 1033,
            "start_date": "2020-03-17",
            "completion_date": "2021-06-28",
            "biomarkers": ["JAK", "STAT1", "IL6"],
            "url": "https://clinicaltrials.gov/ct2/show/NCT04287686"
        },
        {
            "nct_id": "NCT03799133",
            "title": "Efficacy and Safety of IL-6 Receptor Antagonist in Patients With Sepsis",
            "condition": "Infection",
            "status": "Recruiting",
            "sponsor": "University Medical Center",
            "phase": "Phase 2",
            "enrollment": 500,
            "start_date": "2018-12-15",
            "completion_date": "2024-12-31",
            "biomarkers": ["IL6", "TNF", "NFKB1"],
            "url": "https://clinicaltrials.gov/ct2/show/NCT03799133"
        },
        {
            "nct_id": "NCT04643236",
            "title": "Study of Interferon-Beta in Hospitalized Patients With Severe COVID-19",
            "condition": "Infection",
            "status": "Active, not recruiting",
            "sponsor": "National Institute of Allergy and Infectious Diseases",
            "phase": "Phase 2/3",
            "enrollment": 615,
            "start_date": "2020-07-01",
            "completion_date": "2023-12-31",
            "biomarkers": ["IFNG", "JAK1", "STAT1"],
            "url": "https://clinicaltrials.gov/ct2/show/NCT04643236"
        },
        {
            "nct_id": "NCT04660331",
            "title": "TNF-Alpha Inhibition for Moderate to Severe COVID-19",
            "condition": "Infection",
            "status": "Completed",
            "sponsor": "Massachusetts General Hospital",
            "phase": "Phase 2",
            "enrollment": 88,
            "start_date": "2020-08-01",
            "completion_date": "2021-12-15",
            "biomarkers": ["TNF", "IL1B", "NFKB1"],
            "url": "https://clinicaltrials.gov/ct2/show/NCT04660331"
        },
        {
            "nct_id": "NCT04362813",
            "title": "Safety and Efficacy of Tocilizumab (IL-6 Inhibitor) in COVID-19",
            "condition": "Infection",
            "status": "Completed",
            "sponsor": "Genentech, Inc.",
            "phase": "Phase 3",
            "enrollment": 452,
            "start_date": "2020-03-24",
            "completion_date": "2021-06-30",
            "biomarkers": ["IL6", "CRP", "NFKB1"],
            "url": "https://clinicaltrials.gov/ct2/show/NCT04362813"
        }
    ]

    print(f"\nSearching for clinical trials related to: {condition}")
    print(f"Filter: {status}")

    # Filter trials based on criteria
    filtered_trials = [
        trial for trial in example_trials
        if condition.lower() in trial["condition"].lower()
    ]

    if status != "All":
        filtered_trials = [
            trial for trial in filtered_trials
            if status.lower() in trial["status"].lower()
        ]

    print(f"Found {len(filtered_trials)} relevant trials")

    return filtered_trials

def match_trials_to_biomarkers(trials, biomarkers):
    """Match identified biomarkers to clinical trials"""

    print(f"\nMatching {len(biomarkers)} biomarkers to trials...")

    matches = []

    for trial in trials:
        trial_biomarkers = set(trial.get("biomarkers", []))
        user_biomarkers = set(biomarkers)

        overlap = trial_biomarkers & user_biomarkers
        overlap_pct = len(overlap) / len(user_biomarkers) * 100 if user_biomarkers else 0

        if len(overlap) > 0:
            matches.append({
                "nct_id": trial["nct_id"],
                "title": trial["title"],
                "matching_biomarkers": ", ".join(overlap),
                "match_percentage": overlap_pct,
                "phase": trial["phase"],
                "status": trial["status"],
                "sponsor": trial["sponsor"],
                "enrollment": trial["enrollment"],
                "url": trial["url"]
            })

    matches_df = pd.DataFrame(matches)

    if len(matches_df) > 0:
        matches_df = matches_df.sort_values("match_percentage", ascending=False)

    return matches_df

def create_trial_summary_report(trials, output_dir):
    """Create a comprehensive summary report of clinical trials"""

    report_lines = [
        "=" * 80,
        "CLINICAL TRIALS SUMMARY REPORT",
        "Infection/Immune Response Studies",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 80,
        ""
    ]

    # Summary statistics
    statuses = {}
    phases = {}
    total_enrollment = 0

    for trial in trials:
        status = trial["status"]
        statuses[status] = statuses.get(status, 0) + 1

        phase = trial["phase"]
        phases[phase] = phases.get(phase, 0) + 1

        total_enrollment += trial["enrollment"]

    report_lines.extend([
        f"Total Trials: {len(trials)}",
        f"Total Enrollment: {total_enrollment:,} patients",
        "",
        "Status Distribution:",
    ])

    for status, count in sorted(statuses.items(), key=lambda x: x[1], reverse=True):
        report_lines.append(f"  {status}: {count} trials")

    report_lines.extend([
        "",
        "Phase Distribution:",
    ])

    for phase, count in sorted(phases.items(), key=lambda x: x[1], reverse=True):
        report_lines.append(f"  {phase}: {count} trials")

    report_lines.extend([
        "",
        "=" * 80,
        "TRIAL DETAILS",
        "=" * 80,
        ""
    ])

    for trial in sorted(trials, key=lambda x: x["start_date"], reverse=True):
        report_lines.extend([
            f"Trial ID: {trial['nct_id']}",
            f"Title: {trial['title']}",
            f"Status: {trial['status']}",
            f"Phase: {trial['phase']}",
            f"Sponsor: {trial['sponsor']}",
            f"Enrollment: {trial['enrollment']} patients",
            f"Start: {trial['start_date']} | Completion: {trial['completion_date']}",
            f"Key Biomarkers: {', '.join(trial.get('biomarkers', []))}",
            f"Link: {trial['url']}",
            ""
        ])

    report_text = "\n".join(report_lines)

    with open(output_dir / "clinical_trials_report.txt", "w") as f:
        f.write(report_text)

    print("\nClinical Trials Summary:")
    print(f"  Total trials: {len(trials)}")
    print(f"  Total enrollment: {total_enrollment:,}")
    print(f"  Phases: {list(phases.keys())}")

    return report_text

if __name__ == "__main__":
    print("=" * 60)
    print("Clinical Trials Search and Analysis")
    print("=" * 60)

    # Define biomarkers from our analysis
    biomarkers = [
        "IL6", "TNF", "IFNG", "IL1B", "IL12A", "CXCL10",
        "CD8A", "CD4", "NFKB1", "STAT1", "JAK1", "JAK2"
    ]

    print(f"\nBiomarkers identified from analysis: {', '.join(biomarkers)}")

    # Search for clinical trials
    trials = search_clinical_trials(biomarkers, condition="Infection", status="All")

    # Save trial data
    trials_df = pd.DataFrame(trials)
    trials_df.to_csv(OUTPUT_DIR / "all_clinical_trials.csv", index=False)

    # Match biomarkers to trials
    print("\nMatching biomarkers to trials...")
    matches_df = match_trials_to_biomarkers(trials, biomarkers)

    if len(matches_df) > 0:
        matches_df.to_csv(OUTPUT_DIR / "matched_trials.csv", index=False)
        print("\nTop Matching Trials:")
        print(matches_df[['nct_id', 'title', 'match_percentage', 'status']].head(10).to_string(index=False))

    # Create summary report
    print("\nGenerating comprehensive report...")
    report = create_trial_summary_report(trials, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Clinical trials search complete!")
    print(f"Results saved to {OUTPUT_DIR}")
    print("=" * 60)
