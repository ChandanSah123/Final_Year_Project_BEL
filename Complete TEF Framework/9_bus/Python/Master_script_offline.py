import runpy
import traceback
from pathlib import Path


SCRIPT_SEQUENCE = [
    ("Dynamic simulation and data export", "WECCdata1.py"),
    ("Y-bus extraction and Kron reduction", "Y_Pre_F_PF_Krons2.py"),
    ("Critical clearing time search", "Critical_Clearing.py"),
]


def run_stage(stage_name, script_name, base_dir):
    script_path = base_dir / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Required script not found: {script_path}")

    print(f"\n{'=' * 72}")
    print(f"STAGE: {stage_name}")
    print(f"SCRIPT: {script_path.name}")
    print(f"{'=' * 72}")
    runpy.run_path(str(script_path), run_name="__main__")


def main():
    base_dir = Path(__file__).resolve().parent

    for stage_name, script_name in SCRIPT_SEQUENCE:
        run_stage(stage_name, script_name, base_dir)

    print("\nAll offline workflow stages completed successfully.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nMaster workflow failed: {exc}")
        traceback.print_exc()
        raise
