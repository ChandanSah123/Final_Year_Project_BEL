try:
    import psse3603  # type: ignore
    import psspy     # type: ignore
    import dyntools  # type: ignore
    pss_available = True
except Exception:
    pss_available = False
import matplotlib.pyplot as plt
import pandas as pd
import scipy.io
import numpy as np
import os
import re
import config

# ==============================================================================
# 1. SETUP & CONFIGURATION
# ==============================================================================
out_file = config.OUT_FILE

# Initialize PSS/E (safe)
try:
    _i = psspy.getdefaultint()
    _f = psspy.getdefaultreal()
    psspy.psseinit(50)
except Exception:
    _i = 0
    _f = 0.0


# ==============================================================================
# 3. PLOTTING
# ==============================================================================

if not os.path.exists(out_file):
    print(f"CRITICAL ERROR: Output file not found at {out_file}")
else:
    chnf_obj = dyntools.CHNF(out_file)
    short_title, chanid, chandata = chnf_obj.get_data()
    time_data = chandata['time']

    print("Generating Plots...")

    # Plot 1: Rotor Angles (Channels 1, 2, 3)
    plt.figure(figsize=(10, 6))
    for ch in range(1, 5):
        if ch in chandata:
            clean_label = str(chanid[ch]).strip()
            plt.plot(time_data, chandata[ch], label=f"Ch{ch}: {clean_label}")
    plt.title(f"Rotor Angles )")
    plt.ylabel("Angle (Deg)")
    plt.grid(True)
    plt.legend()
    plt.show()


print("--- Script Finished ---")