# file funs.R
# copyright (C) 2022-2025 Artur Gramacki and Jarosław Gramacki
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# ///////////////////////////////////////////////////////////////////////////////////////////
filters_coeff <- function(
    fs = 256, 
    notch = c(49, 51), 
    lowpass = 30, 
    highpass = 1,
    low_cut  = 0.5,
    high_cut = 40) {
  # https://openbci.com/forum/index.php?p=/discussion/193/50hz-notch-filter-coefficients

  ## 50 Hz notch filter
  bf.notch <- butter(2, notch / (fs / 2), "stop")
  freqz(bf.notch)

  # Low pass IIR Butterworth, cutoff at 'lowpass' Hz
  bf.low <- butter(4, lowpass / (fs / 2), "low")
  freqz(bf.low)

  
  # High pass IIR Butterwoth, cutoff at 'highpass' Hz
  bf.high <- butter(4, highpass / (fs / 2), "high")
  freqz(bf.high)
  
  # Bandpass filter IIR Butterworth
  bf.bandpass <- butter(4, c(low_cut, high_cut) / (fs / 2), type = "pass")

  
  list(
    bf.notch = bf.notch, 
    bf.low = bf.low, 
    bf.high = bf.high,
    bf.bandpass = bf.bandpass)
}

# ///////////////////////////////////////////////////////////////////////////////////////////
minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}




