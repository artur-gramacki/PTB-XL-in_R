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
filters.coeff <- function (
    fs = 256,
    notch = c(49, 51), notch.order = 2,
    lowpass = 30, lowpass.order = 4,
    highpass = 1, highpass.order = 4,
    bandpass = c(0.5, 40), bandpass.order = 4,
    bandstop = c(0.5, 40), bandstop.order = 4)
{
  ## Notch filter
  notch <- butter(notch.order, notch / (fs / 2), "stop")
  
  # Low pass IIR Butterworth, cutoff at 'lowpass' Hz
  low <- butter(lowpass.order, lowpass / (fs / 2), "low")
  
  # High pass IIR Butterwoth, cutoff at 'highpass' Hz
  high <- butter(highpass.order, highpass / (fs / 2), "high")
  
  # Bandpass filter IIR Butterworth
  bandpass <- butter(bandpass.order, bandpass / (fs / 2), type = "pass")
  
  # Bandstop filter IIR Butterworth
  bandstop <- butter(bandstop.order, bandstop / (fs / 2), type = "stop")
  
  list(
    notch = notch,
    low = low,
    high = high,
    bandpass = bandpass,
    bandstop = bandstop)
}

# ///////////////////////////////////////////////////////////////////////////////////////////
minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}




