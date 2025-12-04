# automates running TestSIR.exe with different configurations to determine error scaling

import logging
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import time

# COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", ]

logfile = "./model_test.log"
executable_path = "./TestSIR.exe"

with open(logfile, "w"): # clear log file contents
    pass
logging.basicConfig(filename=logfile, level=logging.INFO)
logger = logging.getLogger(__name__)

def run_model(nt: int, nx: int, tf: float):
    logger.info(f"Running with {nt=}, {nx=}, {tf=}")
    start_time = time.perf_counter()
    result = subprocess.run([executable_path, f"--nt={nt}", f"--nx={nx}", f"--tf={tf}", "--test"], capture_output=True)
    end_time = time.perf_counter()
    output = result.stdout.decode()
    logger.info("STDOUT:\n"+output)
    if result.stderr:
        logger.error("STDERR:\n"+result.stderr.decode())
    logger.info(f"Return code: {result.returncode}")
    outstrs = output.split("\n")
    error_norm = None
    for s in outstrs:
        if s.startswith("error norm"):
            if error_norm:
                logger.warning("Multiple err norm candidates")
            try:
                error_norm = float(s.split("=")[1][1:])
            except Exception as e:
                logger.warning(f"Could not parse err norm (string: '{s}', exception: {e})")
    if not error_norm:
        logger.warning("No err norm output found")
        return pd.NA
    logger.info(f"L2 err norm: {error_norm}")
    time_elapsed = end_time - start_time
    logger.info(f"Time to run: {time_elapsed} seconds")
    return error_norm, time_elapsed

def plot_result(df):
    tfs = df.index.unique(level=2)
    nrows = len(tfs)
    fig, axs = plt.subplots(ncols=2, nrows=nrows, sharex=True, sharey=True)

    for i, (tf, _df) in enumerate(df.groupby(level=2)):
        axT = axs[i][0] if nrows > 1 else axs[0]
        axT2 = axT.twinx()
        axX = axs[i][1] if nrows > 1 else axs[1]
        axX2 = axX.twinx()

        axT.annotate(r"$T_f={}$".format(tf),
            xy=(0, 0.5),
            xytext=(-axT.yaxis.labelpad - 5, 0),
            xycoords=axT.yaxis.label,
            textcoords='offset points',
            size='large',
            ha='right',
            va='center')

        for nx, _dfT in _df.groupby(level=1):
            result_nts = []
            result_errs = []
            result_times = []
            for (nt, _, _), v in _dfT.iterrows():
                result_nts.append(nt)
                result_errs.append(v["err"])
                result_times.append(v["time"])
            axT.plot(result_nts, result_errs, label=f"{nx=}", marker="o")
            axT2.plot(result_nts, result_times, label=f"t, {nx=}", marker="o", linestyle="dashed")

        for nt, _dfX in _df.groupby(level=0):
            result_nxs = []
            result_errs = []
            result_times = []
            for (_, nx, _), v in _dfX.iterrows():
                result_nxs.append(nx)
                result_errs.append(v["err"])
                result_times.append(v["time"])
            axX.plot(result_nxs, result_errs, label=f"{nt=}", marker="o")
            axX2.plot(result_nxs, result_times, label=f"t, {nt=}", marker="o", linestyle="dashed")

        axT.legend()
        # axT2.legend()
        # axT.set_xscale("log", base=2)
        axT.set_xlabel("nt")
        # axT.set_yscale("log", base=10)
        axT.set_ylabel("L2 Error Norm")
        
        axX.legend()
        # axX2.legend()
        # axX.set_xscale("log", base=2)
        axX.set_xlabel("nx")
        # axX.set_yscale("log", base=10)
        axX2.set_ylabel("Time Taken (s)")

    fig.tight_layout()
    plt.show()

def save_result(df, filepath):
    df.to_csv(filepath)

def test_scaling(nts, nxs, tfs):
    index = pd.MultiIndex.from_product([nts, nxs, tfs], names=["nt","nx","tf"])
    df = pd.DataFrame(data=pd.NA, columns=["err", "time"], index=index)
    for i in df.index:
        err, time = run_model(*i)
        df.loc[i, "err"] = err
        df.loc[i, "time"] = time
    return df

if __name__ == "__main__":
    # TFS = [1.0]
    # NUMS = [8, 16]
    NUMS = [16, 24, 32, 64, 128]
    TFS = [1.0, 2.0]
    result = test_scaling(NUMS, NUMS, TFS)
    # result = pd.DataFrame(data=[[1+i, 18.5-i] for i in range(18)], columns=["err", "time"], index=pd.MultiIndex.from_product([[16, 24, 32], [16, 24, 32], [1.0, 2.0]], names=["nt","nx","tf"]))
    save_result(result, "test_scaling_results.csv")
    plot_result(result)
    
