"""run_from_config.py
Run any project script from one place using config.WORK_DIR as cwd.
Usage examples:
    python run_from_config.py --script WSCCPLOT
    python run_from_config.py --script Y_Pre_F_PF_Krons2 --args "--some-arg val"
    python run_from_config.py --list
"""
import argparse
import subprocess
import sys
import os
import config

SCRIPTS = {
    'WSCCPLOT': 'WSCCPLOT.py',
    'WSCCDATA1': 'WSCCdata1.py',
    'CRITICAL': 'Critical_Clearing.py',
    'COMBINED_CCT': 'combined_CCT.py',
    'Y_PRE_KRONS2': 'Y_Pre_F_PF_Krons2.py',
}


def run_script(script_key, extra_args=None):
    if script_key not in SCRIPTS:
        raise SystemExit(f"Unknown script '{script_key}'. Use --list to see valid names.")

    script_path = os.path.join(config.WORK_DIR, SCRIPTS[script_key])
    if not os.path.exists(script_path):
        raise SystemExit(f"Script not found: {script_path}")

    cmd = [sys.executable, script_path]
    if extra_args:
        # allow passing additional args as a single string
        cmd += extra_args.split()

    print(f"Running: {' '.join(cmd)} (cwd={config.WORK_DIR})")
    proc = subprocess.run(cmd, cwd=config.WORK_DIR)
    return proc.returncode


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run project scripts from one place (uses config.WORK_DIR).')
    parser.add_argument('--list', action='store_true', help='List available scripts')
    parser.add_argument('--script', type=str, help='Script key to run (case-insensitive)')
    parser.add_argument('--args', type=str, default='', help='Extra arguments to pass to the script')

    args = parser.parse_args()

    if args.list:
        print('Available scripts:')
        for k, v in SCRIPTS.items():
            print(f" - {k}: {v}")
        sys.exit(0)

    if not args.script:
        parser.print_help()
        sys.exit(1)

    key = args.script.strip().upper()
    rc = run_script(key, args.args)
    sys.exit(rc)
