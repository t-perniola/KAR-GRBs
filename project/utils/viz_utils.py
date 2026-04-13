import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.tsa.stattools import acf

# -------------------------------------------------------
# VISUALIZATION UTILITIES
# -------------------------------------------------------
def viz_df_cols(satellite_df, cols_to_plot, n_rows=4, n_cols=2):
  # Get the list of columns to plot
  n_cols_to_plot = len(cols_to_plot)
  fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(12, n_rows * 2.5)) # Adjust figure height based on nrows
  for i, feature_col in enumerate(cols_to_plot):
    row, col = divmod(i, n_cols)  # Compute subplot position
    satellite_df.plot.scatter(
          x='timestamp', y=feature_col, s=2, title=f'{feature_col} over time', ax=ax[row, col]
      )
  # Hide any unused subplots
  for i in range(n_cols_to_plot, n_rows * n_cols):
      row, col = divmod(i, n_cols)
      fig.delaxes(ax[row, col])
  plt.tight_layout()
  plt.show()

def viz_lightcurve(satellite_df, time_res, det=1, interval = (0, 150)):
  start, end = interval
  plt.figure(figsize = (12, 5))
  plt.step(x = satellite_df['timestamp'][start:end], y = satellite_df[f'rates_det_{det}'][start:end])
  plt.xlabel("Time")
  plt.ylabel(f"Photon count rate (detector {det})")
  plt.title(f"{time_res}-second binned light curve")
  plt.grid(alpha=0.3)
  plt.show()

def viz_countrates(
    satellite_df,
    rates,
    time_res,
    count_str,
    title,
    y_zoom=None,
):
    detector_names = [
        "Z1",
        "Z0",
        "X1",
        "X0",
        "Y1",
        "Y0",
    ]

    plt.figure(figsize=(12, 7))
    for i in range(len(rates[0])):
        plt.scatter(
            x=satellite_df["timestamp"],
            y=satellite_df[f"{count_str}{i+1}"],
            s=1,
            label=detector_names[i],
        )

    plt.title(title)
    plt.xlabel(f"Time (bins of {time_res}s)")
    plt.ylabel("Counts")

    if y_zoom is not None:
        plt.ylim(*y_zoom)
    else: plt.ylim(bottom=1000)

    plt.legend(markerscale=4)
    plt.grid(alpha=0.3)
    plt.show()



def viz_df_heatmap(satellite_df, count_str, time_res="15"):
  # Select only the detectors' count rate columns
  detector_cols = [col for col in satellite_df.columns if f"{count_str}" in col]
  counts = satellite_df[detector_cols]
  # Standardize each detector column: (x - mean) / std
  # So values below mean are negative (blue), above mean positive (red)
  counts_normalized = (counts - counts.mean()) / counts.std()
  # Plot heatmap
  plt.figure(figsize=(12, 5))
  sns.heatmap(
      counts_normalized.T,  # transpose → detectors on y-axis, time on x-axis
      cmap=sns.color_palette("RdBu_r", as_cmap=True),   # red = above mean, blue = below mean, white = near mean
      center=0,             # zero centered on white
      cbar_kws={'label': 'Deviation from Mean (z-score)'}
  )
  plt.title(f"Detector counts over time — time_res: {time_res}")
  plt.xlabel("Time")
  plt.ylabel("Detectors")
  plt.tight_layout()
  plt.show()

def viz_histograms(satellite_df, bins=20):
  nrows = 3
  fig, ax = plt.subplots(nrows=nrows, ncols=nrows, figsize=(12, 7))
  fig.suptitle("Histogram of columns")
  # Get the list of columns to plot
  cols_to_plot = satellite_df.columns[1:10]
  for i, col in zip(range(len(cols_to_plot)), cols_to_plot):
      ax[i // nrows, i % nrows].hist(satellite_df[col], bins=bins)
      ax[i // nrows, i % nrows].set_title(col)
      ax[i // nrows, i % nrows].set_xlabel(col)
      ax[i // nrows, i % nrows].set_ylabel("Frequency")
  plt.tight_layout()
  plt.show()

def viz_kdes(satellite_df, cols_to_plot, n_rows=4, n_cols=2):
  fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(18, 10))
  fig.suptitle("KDE plots of columns")
  # Get the list of columns to plot (first 8 columns from index 1)
  for i, feature in enumerate(cols_to_plot):
      row, col = divmod(i, n_cols)
      sns.kdeplot(data=satellite_df[feature], ax=ax[row, col])
      ax[row, col].set_title(feature)
  plt.tight_layout()
  plt.show()

def viz_ac_plots(satellite_df, cols_to_plot, n_rows=4, n_cols=2):
  # AC plot for input features
  fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(18, 10))
  fig.suptitle("Autocorrelation plots of the input features")
  # Get the list of columns to plot (first 8 columns from index 1)
  for i, feature in enumerate(cols_to_plot):
      row, col = divmod(i, n_cols)
      acf_values = acf(satellite_df[feature], nlags=3000)
      ax[row, col].plot(acf_values)
      ax[row, col].set_title(feature)
  plt.tight_layout()
  plt.show()
