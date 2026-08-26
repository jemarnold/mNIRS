#' @name moxy_ramp.xlsx
#'
#' @title 2 Hz PerfPro export of Moxy data
#'
#' @description Exported from PerfPro Studio software, recorded
#'   at 0.5 Hz no smoothing and exported at 2 Hz. Containing a
#'   ramp incremental cycling protocol, placed on bilateral
#'   vastus lateralis muscle sites. Intentional data errors
#'   (outliers, invalid values, and missing `NA` values) have
#'   been introduced to demonstrate *{mnirs}* cleaning
#'   functions.
#'
#' @docType data
#'
#' @format .xlsx file with five columns and 2202 rows:
#'   \describe{
#'     \item{mm-dd}{Recording date (dd-MMM format).}
#'     \item{hh:mm:ss}{Time of day (hh:mm:ss format).}
#'     \item{SmO2 Live}{Muscle oxygen saturation, left leg (%). Contains
#'       simulated erroneous and missing samples.}
#'     \item{SmO2 Live(2)}{Muscle oxygen saturation, right leg (%).}
#'     \item{Lap}{Lap marker (integer).}
#'   }
#'
#'   Channel mapping for [read_mnirs()]:
#'   - `nirs_channels = c("SmO2 Live", "SmO2 Live(2)")`
#'   - `time_channel = c("hh:mm:ss")`
#'   - `event_channel = c("Lap")`
#'   - `interval_times = list(start = c(204, 878))` (start and end of exercise)
#'
#' @source Moxy Monitor (Fortiori Design LLC), exported via PerfPro Studio
#'   desktop software (https://perfprostudio.com/).
#'
#' @seealso [read_mnirs()], [example_mnirs()]
#'
#' @examples
#' example_mnirs("moxy_ramp")
#'
NULL


#' @name moxy_intervals.csv
#'
#' @title 0.5 Hz Moxy onboard export
#'
#' @description Exported from Moxy onboard recording at 0.5 Hz
#'   no smoothing. Containing four 4-minute cycling work intervals,
#'   placed on the vastus lateralis muscle site.
#'
#' @docType data
#'
#' @format .csv file with seven columns and 936 rows:
#'   \describe{
#'     \item{mm-dd}{Recording date (dd-MMM format).}
#'     \item{hh:mm:ss}{Recording time of day (hh:mm:ss format).}
#'     \item{SmO2 Live}{Muscle oxygen saturation, raw signal (%).}
#'     \item{SmO2 Averaged}{Muscle oxygen saturation, rolling average (%).}
#'     \item{THb}{Total haemoglobin (arbitrary units).}
#'     \item{Lap}{Lap marker (integer). Not typically in use.}
#'     \item{Session Ct}{Session count of recordings.}
#'   }
#'
#'   Channel mapping for [read_mnirs()]:
#'   - `nirs_channels = c("SmO2 Live", "SmO2 Averaged", "THb")`
#'   - `time_channel = c("hh:mm:ss")`
#'   - `interval_times = list(start = c(124, 486, 848, 1210), end = c(364, 726, 1088, 1450))`
#'
#' @source Moxy Monitor (Fortiori Design LLC), exported via Moxy Portal App.
#'   (https://www.moxymonitor.com/)
#'
#' @seealso [read_mnirs()], [example_mnirs()]
#'
#' @examples
#' example_mnirs("moxy_intervals")
#'
NULL


#' @name train.red_intervals.csv
#'
#' @title 10 Hz Train.Red App export
#'
#' @description Exported from Train.Red app, recorded at 10 Hz.
#'   Containing two 5-minute cycling work intervals, placed on
#'   bilateral vastus lateralis muscle sites. Some data channels
#'   have been omitted to reduce file size.
#'
#' @docType data
#'
#' @format .csv file with header metadata and 10 columns and 11995 rows:
#'   \describe{
#'     \item{Timestamp (seconds passed)}{Elapsed time (s).}
#'     \item{Lap/Event}{Lap number (numeric).}
#'     \item{SmO2}{Muscle oxygen saturation, filtered (%). Two channels have
#'       duplicated names. If both are called, the second will be renamed to
#'       `SmO2_1`.}
#'     \item{SmO2 unfiltered}{Muscle oxygen saturation, raw signal (%). Two
#'       channels have duplicated names. If both are called, the second will
#'       be renamed to `SmO2 unfiltered_1`.}
#'     \item{O2HB unfiltered}{Oxyhaemoglobin concentration, raw signal
#'       (arbitrary units). Two channels have duplicated names. If both are
#'       called, the second will be renamed to `O2HB unfiltered_1`.}
#'     \item{HHB unfiltered}{Deoxyhaemoglobin concentration, raw signal
#'       (arbitrary units). Two channels have duplicated names. If both are
#'       called, the second will be renamed to `HHb unfiltered_1`.}
#'   }
#'
#'   Channel mapping for [read_mnirs()]:
#'   - `nirs_channels = c("SmO2", "SmO2 unfiltered", "O2HB unfiltered", "HHb unfiltered")`
#'   - `time_channel = c("Timestamp (seconds passed)")`
#'   - `event_channel = c("Lap/Event")`
#'   - `interval_times = list(start = c(2150.09, 2872.28), end = c(2452.26, 3167.98))`
#'   - `interval_times = list(start = c(65.94, 788.13), end = c(368.11, 1083.83))` from zero_time
#'
#' @source Train.Red (Train.Red B.V.), exported via Train.Red app
#'   (https://train.red/)
#'
#' @seealso [read_mnirs()], [example_mnirs()]
#'
#' @examples
#' example_mnirs("train.red")
#'
NULL


#' @name artinis_intervals.xlsx
#'
#' @title 10 Hz Artinis Oxysoft export recorded with Oxymon MKIII
#'
#' @description Exported from Artinis Oxysoft, recorded on Oxymon MKIII at
#'   50 Hz and exported at 10 Hz. Containing two 5-minute cycling work
#'   intervals and an ischaemic occlusion, placed on the vastus lateralis
#'   muscle site.
#'
#' @docType data
#'
#' @format .xlsx file with header metadata and five columns and 20919 rows:
#'   \describe{
#'     \item{Column 1}{Sample index (divide by sample rate for seconds).}
#'     \item{Column 2}{O2Hb: oxyhaemoglobin concentration change (\eqn{\mu}M).}
#'     \item{Column 3}{HHb: deoxyhaemoglobin concentration change (\eqn{\mu}M).}
#'     \item{Column 4}{Event marker (character).}
#'     \item{Column 5}{Unmarked event label (character).}
#'   }
#'
#'   Channels are detected automatically from the file legend, or can be
#'   specified explicitly for [read_mnirs()]:
#'   - `nirs_channels = c(O2Hb = 2, HHb = 3)`
#'   - `time_channel = c(sample = 1)`
#'   - `event_channel = c(event = 4)`
#'   - `interval_times = list(start = c(158, 999, 1750) end = c(493, 1333, 1961))` two intervals and post-exercise occlusion
#'
#' @source Artinis Medical Systems. Oxymon MKIII, exported via Oxysoft desktop
#'   software (https://artinis.com/)
#'
#' @seealso [read_mnirs()], [example_mnirs()]
#'
#' @examples
#' example_mnirs("artinis_intervals")
#'
NULL


#' @name portamon-oxcap.xlsx
#'
#' @title 10 Hz Artinis Oxysoft export recorded with Portamon
#'
#' @description Exported from Artinis Oxysoft, recorded on
#'   Portamon at 10 Hz on the vastus lateralis muscle. Containing
#'   two trials of repeated occlusion oxidative capacity testing,
#'   each with 17 occlusions.
#'
#' @docType data
#'
#' @format .xlsx file with header metadata and six columns and 7943 rows:
#'   \describe{
#'     \item{Column 1}{Sample index (divide by sample rate for seconds).}
#'     \item{Column 2}{tHb: total haemoglobin concentration change
#'       (\eqn{\mu}M).}
#'     \item{Column 3}{HHb: deoxyhaemoglobin concentration change (\eqn{\mu}M).}
#'     \item{Column 4}{O2Hb: oxyhaemoglobin concentration change (\eqn{\mu}M).}
#'     \item{Column 5}{Event marker (character).}
#'     \item{Column 6}{Unmarked event label (character).}
#'   }
#'
#'   Channels are detected automatically from the file legend (the unmarked
#'   label column is named `"labels"`), or can be specified explicitly for
#'   [read_mnirs()]:
#'   - `nirs_channels = c(THb = 2, HHb = 3, O2Hb = 4)`
#'   - `time_channel = c(sample = 1)`
#'   - `event_channel = c(event = 5, label = "col_6")`
#'
#' @source Artinis Medical Systems. Portamon, exported via Oxysoft desktop
#'   software (https://artinis.com/)
#'
#' @seealso [read_mnirs()], [example_mnirs()]
#'
#' @examples
#' example_mnirs("portamon")
#'
NULL
