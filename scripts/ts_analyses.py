#!/usr/bin/python
#
# Tools to work out climate time series
# Creation:       February 2017
# Contributor(s): F. Massonnet, UCL
def eftimescale(data, \
                period      = None, \
                printfig    = None)  :

  import sys
  import numpy as np
  
  """ Input:  - data : 1-D numpy array of size N containing a time series. 
                The series is assumed to be sampled at regular intervals (even 
                time spacing). For meaningful analyses the data should have
                been detrended appropriately. Read the "Important" note below.

              - period: integer giving a possible periodicity in the persistence
                        of the data. If there is no a priori reason to believe
                        that this is the case, specify None or nothing.
 
              - printfig is a boolean to print a summary figure of the auto
                correlation function

      Outputs: 1. "tau", a Numpy array of size 1 or size "period", depending on whether
                   a period was specified, giving the e-folding time scale of the
                   process. The e-folding time scale is defined as the time it takes
                   for the auto-correlation of the process to drop below 1 / e.
                   Linear interpolation between the successive lags for which the empirical auto-
                   correlation function is above, and below 1/e is conducted.

               2. "acf". If period = None, acf is a numpy array giving the ACF until
                  the lag is too long to give meaningful correlations.
                         If period is not None, a list of numpy arrays (one for each
                         element of the period) 
      
      Method: - The empirical auto-correlation function (acf) is first
                estimated by correlating the signal with itself at various
                lags. Then, the time where acf is less than 1 / e for the first
                time is determined.
                
      Important: 1. This fit assumes an AR1 structure, i.e. it returns the time
                    scale that would make an AR1 process have the same auto-
                    correlation function as the input data.

                 2. The e-folding time scale characterizes the "fast", "short"
                    memory of the process as it indicates the rate at which
                    the auto-correlation drops towards zero, but tells nothing
                    about what could happen at longer lags. It shouldn't be over-
                    interpreted!

                 3. Nearly all climate time series are composed of a forced 
                    component and a a random component (noise, possibly with
                    structure). Usually one attempts to remove the forced
                    component from the raw signal prior to estimation of 
                    persistence, because one is interested in the intrinsic
                    serial memory of the time series and not to that associated
                    with the forcing (that is known). This task is left at the 
                    discretion of the user: the current function does not
                    perform any detrending. A function "detrend" in this
                    file can do that.
                 
      Author: F. Massonnet (UCL-BSC)
      Date: Feb 22., 2017
            June 2017: update for periodicity
  """

  # Basic checks: does the data look as expected?
  if len(data.shape) != 1:
    sys.exit("(eftimescale): non-conform input data")

  N = len(data)
  if N <= 2:
    sys.exit("(eftimescale): input data too short")

  if sum(1.0 * np.isnan(data)) > 0:
    # If there is nan, exit
    sys.exit("(eftimescale): input data contains nan(s)!")

  if period is None:
    #-------------------------------
    # Case 1: No period is specified
    #-------------------------------
    #
    # We correlate the signal with itself at all lags, starting from 0, and 
    # until we reach a value less than 1 / e 

    acf = list()     # We will add item per item
    reached_threshold = False # Will be false until correlation < 1/e
    keep_going = True # To continue compute ACF until sample size is too small
    lag = 0          # Will be incremented

    while keep_going:

      predictand  = data[lag:]
      predictor   = data[:N - lag]
      
      # The data may contain NaNs. Delete them at the same place in both series
      X = predictor
      Y = predictand
 
      if len(X) <= 3 or len(Y) <= 3:
        # If the lag is too large, then obviously the number of data in the
        # predictor and predictand will become small. In that case, abort.
        break 

      # Compute correlation between predictor and predictand
      corr = np.corrcoef(X, Y)[0][1]
 
      # If, for the first time, correlation drops below 1/e, then record the
      # e-folding time scale and say that we reached that point
      if not reached_threshold and corr < 1.0 / np.exp(1):
        reached_threshold = True
        # Interpolate between previous and current lag to determine the exact crossing of 1/e
        tau = lag  + ((lag - 1.0) - lag) * (1.0 / np.exp(1) - corr) / (acf[-1] - corr)
      
      # Append the correlation to the auto-correlation function
      acf.append(corr)
  
      # Increment the lag
      lag += 1

    acf = np.asarray(acf) # Reformat as a numpy array

  else:
    # -----------------------------
    # Case 2: A period is specified
    # -----------------------------
    #
    # For explanation of how we do, assume a 4-periodic signal (e.g. seasons)
    # running over 4 years
    # data = [Wi0 Sp0 Su0 Fa0 Wi1 Sp1 Su1 Fa1 Wi2 Sp2 Su2 Fa2 Wi3 Sp3 Su3 Fa3 Wi4 Sp4 Su4 Fa4]
    # (Wi:Wintner; Sp: Spring; Su: Summer; Fa: Fall)
    # We assume the data is already detrended
 
    # The idea is to iterate for each period element: Wi, Sp, Su, Fa i.e. indices 0, 1, 2, 3
    # and correlate all indices with the next one, then to above, etc:
    # 1. Winter (jp = 0)
    # ------------------
    #   at lag 0:
    #     corr([Wi0, Wi1, Wi2, Wi3, Wi4], [Wi0, Wi1, Wi2, Wi3, Wi4])
    #   at lag 1:
    #     corr([Wi0, Wi1, Wi2, Wi3, Wi4], [Sp0, Sp1, Sp2, Sp3, Sp4])
    #   at lag 2:
    #     corr([Wi0, Wi1, Wi2, Wi3, Wi4], [Su0, Su1, Su2, Su3, Su4])
    #   ...
    #   until the correlation is less than 1/e. We record the time scale
    #   for that period element jp = 0
    # 2. Spring (jp = 1)
    # ------------------
    #   at lag 0:
    #     corr([Sp0, Sp1, Sp2, Sp3, Sp4], [Sp0, Sp1, Sp2, Sp3, Sp4])
    #   etc...

    acf = list()
    
    tau = np.empty([period])
    tau[:] = np.nan

    for jp in range(period):
      # For each period
      acf.append(list()) # We create a list that will contain all correlations at different lags
                         # This list will eventually be converted to a numpy array

      print("(eftimescale) Computing auto-correlations " + str(jp + 1) + " / " + str(period))

      predictor = data[np.arange(jp, N, period)] # e.g. all winters

      lag = 0
      keep_going = True         # to maintain the loop over lags
      reached_threshold = False # whether the current autocorrelation is < 1/e
 
      while keep_going:

        predictand = data[np.arange(jp + lag, N, period)]

        # After one period the predictand will become shorter than the predictor
        # because of the lag.
        # For our seasonal example, at lag 5:
        # Predictor =  [Wi0, Wi1, Wi2, Wi3, Wi4] 
        # Predictand = [Sp1, Sp2, Sp3, Sp4]

        # Make sure the predictor and predictand have same
        # length

        predictor = predictor[range(len(predictand))]

        X = predictor
        Y = predictand

        if len(X) <= 3 or len(Y) <=3: # if lag is too large,
                                      # not enough data to compute correlation
          break                       # --> exit the loop

        corr = np.corrcoef(X, Y)[0][1]

        # If reached less than 1/e for the first time
        if not reached_threshold and corr < 1.0 / np.exp(1):
          reached_threshold = True
          # Interpolate between previous and current lag to determine the exact crossing of 1/e
          tau[jp] = lag  + ((lag - 1.0) - lag) * (1.0 / np.exp(1) - corr) / (acf[jp][-1] - corr)

        # Append the result
        acf[jp].append(corr)
  
        # Increment the lag
        lag += 1

      # We have gone out of the loop because lag was too long.
      # Convert the data to a numpy array, and we are done
      acf[jp] = np.asarray(acf[jp])


  if printfig is not None:
    import matplotlib.pyplot as plt
    plt.close("all")
    if period is None:
      # Print a single figure showing the signal and the decorrelation
      fig = plt.figure("show", figsize = (8, 8))

      plt.subplot(211)
      time = np.arange(len(data))
      plt.plot(time, data, 'k.')
      plt.xlabel("Time")
      plt.ylabel("signal")

      plt.subplot(212)
      time = np.arange(len(acf))
      plt.plot(acf,  'b.')
      plt.xlabel("Lag")
      plt.ylabel("r")
    
      plt.plot([0, 1e9], [1.0 / np.exp(1), 1.0 / np.exp(1)], 'k--')
      plt.plot([tau, tau], [-1e9, 1e9], 'k--')
      plt.plot(tau, 0, 'r.', ms = 10)

      plt.plot([0, 1e9], [0, 0], 'k')
      plt.xlim(0, 6 * tau)
      plt.ylim(-0.5, 1.0)

      plt.tight_layout()
      plt.savefig(printfig)
      print("(eftimescale): " + printfig + " printed")
      plt.close("all")
  
    else:
      # Print a 4 by 3 figure where 12 initial  times have been sampled through
      # the period
      fig = plt.figure("acf", figsize = (10, 10))
      jS = 1 # for subplot
      indices = np.linspace(0, period - 1, 12)
      for jp in [np.int(i) for i in indices]:
        time = np.arange(len(acf[jp, :]))

        fig.add_subplot(3, 4, jS)
        plt.title("Predictor:\nperiod index " + str(jp))

        plt.plot(time, acf[jp, :], 'b')
        plt.plot(time, acf[jp, :], 'b.')
        plt.xlabel("Lag")
        plt.ylabel("r")

        plt.plot([0, 1e9], [1.0 / np.exp(1), 1.0 / np.exp(1)], 'k--')
        plt.plot([tau[jp], tau[jp]], [-1e9, 1e9], 'k--')

        plt.plot(tau[jp], 0, 'r.', ms = 10)

        plt.plot([0, 1e9], [0, 0], 'k-')
        plt.xlim(0, 6 * np.max(tau))
        plt.ylim(-0.5, 1)
        jS += 1
      plt.tight_layout()
      figname = "acf.png"
      plt.savefig(figname)
      print("(eftimescale): " + figname + " printed")
      plt.close("acf")

  return tau, acf

def unleap(dates, data):
  """ 
      Input:   dates  (datetime object) and data (numpy array)
      Output:  dates and data without the elements with such a position
               that dates[element] is a 29th of February
  """

  import numpy as np
  new = [(dates[i], data[i]) for i in range(len(dates)) if (dates[i].month != 2 or dates[i].day != 29)]
  return [i[0] for i in new] , np.asarray([i[1] for i in new])


def detrend(data, order = 1, period = None):
  import numpy as np
  import sys

  """ Input: data: 1-D numpy array of size N, assumed to be sampled at
                   evenly spaced times
             order: order of the polynomial for detrending
             period: possible existing periodicity of the signal, coming e.g.
                     from external forcing,  expressed in units time steps. 
                     That is, data[i] and data[i + period] correspond
                     to two realizations of the process at times where the
                     forcing might be similar. Common examples include
                     the seasonal cycle forcing, or the diurnal forcing.
 
                     If "period" is not None, the detrending is performed
                     separately for each time step (e.g., all 1st of January,
                     all 2nd of January, ..., all 31st of December in case
                     of annual cycle.

                     If "period" is None, the detrending is performed on
                     the given time series
                                          
      Output: the signal detrended using a least-square polynomial
              regression of order "order"

  """ 

  if len(data.shape) != 1:
    sys.exit("(detrend): non-conform input data")

  # Remove possible nans from the data. All the regression (ie polyfit)
  # parameters will be estimated based on the no-nan data but the 
  # residuals will be computed from the original data in order 
  # to keep the same size and to restitute NaNs where they appeared

  data_nonan = data[~np.isnan(data)]

  N       = len(data)
  N_nonan = len(data_nonan)

  # If the signal has no periodicity, we just make a linear regression
  if period is None:
    time_nonan = np.arange(N_nonan)
    time       = np.arange(N)
    p = np.polyfit(time_nonan, data_nonan, order)
    residuals = data - np.sum([p[i] * time ** (order - i) \
                  for i in range(order + 1)], axis = 0)
  
  # If the signal contains a periodical component, we do the regression
  # time step per time step
  else:
    residuals = np.empty([N])  

    # For each time step of the period, detrend
    for jP in np.arange(period):
      raw         = data[np.arange(jP, N, period)]
      raw_nonan   = raw[~np.isnan(raw)]
      time        = np.arange(len(raw))
      time_nonan  = np.arange(len(raw_nonan))
      p = np.polyfit(time_nonan, raw_nonan, order)
      residuals[np.arange(jP, N, period)] = \
            raw - np.sum([p[i] * time ** (order - i) \
                            for i in range(order + 1)], axis = 0)

    # Note that another common option is to first remove a seasonal cycle
    # and then detrend the anomalies. However this assumes that a cycle can
    # be estimated, which in presence of a trend is tricky because the
    # trend component interferes with the mean. I have tried that and it
    # gives ugly step-wise anomalies. Detrending day per day seems the
    # the most natural way to do, at least as long as we assume that the 
    # raw signal at some time is the result of a seasonal cycle depending
    # on the position of the time step in the period, plus a common trend,
    # plus some noise.
    
  return residuals


def polychange(data, time = None, order = 1, printfig = None):
  """ Input:  data, a 1-D numpy array
              time, a 1-D numpy array of same size telling
                    what are the sampling times (can be uneven)
                    If none, assumes evenly spaced times

              order, the degree for the polynomial fit
 
              printfig, a string with the name of the figure to print
                        (series + fit + 1 sigma)

      Output: change, a scalar giving the change from polynomial fitting
                      evaluated at the last minus first time
              sd    , a scalar giving the standard deviation on this 
                      change
  """
  import numpy as np
  import sys

  if time is None:
    time = np.arange(len(data))
    
  if len(data.shape) != 1:
    sys.exit("(polychange) Input data is not 1-D")
 
  nt = len(time)

  # Here a little bit of explanation is necessary.
  # Polynomial regression of order "q" takes the form
  #      T' * A = Y_hat 
  # where
  #          t_1^q  t_1^(q-1) ... t_1  1
  #      T=  t_2^q  t_2^(q-1) ... t_2  1
  #          ...
  #          t_M^q  t_M^(q-1) ... t_M  1
  # and
  #      A = [a_0 a_1 a_2 ... a_q]'
  #     is the vector of coefficients ofthe regression (the ' denotes transpose)
  # and Y_hat is the M-dimensional vector that we hope is close to the data we fit

  # np.polyfit returns "p" and "C", such that
  # p = the A that minimizes the error
  # C   is the covariance matrix of p

  p, C = np.polyfit(time[~np.isnan(data)], data[~np.isnan(data)], order, cov = True)

  # Now, we know through statistics that if A ~ N(mu, Sigma)
  # then b + X * A ~ N(b + X * mu, X * Sigma * X')

  # Thus, Y_hat (the vector resulting from the fit at the sample points) follows
  #
  # Y_hat ~ N( T' * p, T' * C * T)

  # Compute regression (the term T' * p)
  tt = np.linspace(np.nanmin(time), np.nanmax(time), 1000)

  reg = np.sum([p[i] * tt ** (order - i) for i in range(order + 1)], axis = 0)

  # Compute co-variance matrix of Y_hat
  TT = [tt ** (order - i) for i in range(order + 1)]
  covfit = np.dot(np.dot(np.transpose(TT), C), TT)
  onesigma = [np.sqrt(covfit[i, i]) for i in range(len(tt))]

  # Now let's call Y the fit evaluated at the last time and X at the first time.
  # The change associated to the fit is Y - X
  # Its expectation is E(Y-X) = E(Y) - E(X) = differences in the regression
  change_mean = reg[-1] - reg[0]

  # Its variance is V(Y - X) = V(Y) + V(X) - 2 * C(X, Y). 
  # From there, easy to compute the standard deviation
  change_std = np.sqrt(covfit[-1, -1] + covfit[0, 0] - 2 * covfit[0, -1])

  if printfig is not None:
    import matplotlib.pyplot as plt

    fig = plt.figure("polychange")
    plt.plot(time, data, "b.", ms = 10)
    plt.plot(tt, reg, "g", lw = 3)
    plt.plot(tt, [x + 1.96 * y for x, y in zip(reg, onesigma)], "g--", lw = 1)
    plt.plot(tt, [x - 1.96 * y for x, y in zip(reg, onesigma)], "g--", lw = 1)
    plt.savefig(printfig)
    print("(polychange) " + printfig + " printed")
    plt.close("polychange")
  
  return change_mean, change_std

def example(): 
  # In this example we pretend to analyse a time series of Arctic sea ice 
  # volume that we build ourselves.
 
  # ----------------------------------------------------------------
  # NOTE: a seed is used to make the result reproducible, but it can 
  #       of course be switched off
  # ----------------------------------------------------------------
  import numpy as np
  import matplotlib.pyplot as plt
  import datetime
  
  np.random.seed(123)

  plt.close("all")

  a = np.exp(-1.0 / 1000)    # Parameter of AR process:
                            # noise(t+1) = a * noise(t) + e(t)
                            # where e(t) is a gaussian    
  std = 0.1              # with zero mean and sd "std". 
                            # Such a process has an auto-correlation function
                            # f(t) = a^t where t is the time. 
                            # Now a^t = exp(ln(a) * t) = exp(t / (-1 / ln(a)))
                            # So the corresponding time scale is -1 / ln(a)
                            # Put the other way around, specify a value
                            # a = exp(-1 / tau) for a decorrelation time scale
                            # of tau
                            
  period = 365              # periodicity
  N = period * 38           # Number of time steps
  base = datetime.datetime(1978, 12, 31)
  dates = [base + datetime.timedelta(days = int(i)) for i in range(N)]

  cycle = np.empty([N])     # seasonal cycle 

  b = -15.0 / N             # slope for forced trend (units volume
                            # per unit time). Let's say that we loose
                            # 15 thousands km3 over the period

  noise = np.empty([N])
  noise[0] = std * np.random.randn()

  # Building a seasonal cycle first
  for jP in np.arange(period):
    cycle[np.arange(jP, N, period)] = \
            30.0 + 10 * np.sin(2 * np.pi / period * jP)

  # Building a trend component
  trend = b * np.arange(N)
  
  # Building a noise component with AR1 structure
  for jN in np.arange(N - 1):
    noise[jN + 1] = a * noise[jN] + std * np.random.randn()
    
  # Creating the full signal
  X = cycle + trend + noise
  
  # --------------------------------------------------------------
  # Now let's pretend that we are given X, and we want to study it
  # --------------------------------------------------------------
  # We know that X has a periodic component and a linear trend
  X_d  = detrend(X, period = period)

  # The time scale is estimated on the detrended signal, of course
  tau, acf = eftimescale(X_d, printfig = "synthetic.png")

  tau_noise, acf_noise = eftimescale(noise) 
                                 # The time scale estimated from
                                 # the noise component

  tau_true = -1 / np.log(a)      # The true time scale

  print("Time scale estimated from the signal:".ljust(40) + \
           str(np.round(tau, 2)).ljust(8) + " units time steps")
  print("Time scale if noise term was known:".ljust(40)   + \
           str(np.round(tau_noise, 2)).ljust(8) + " units time steps")
  print("True time scale:".ljust(40) + \
           str(np.round(tau_true, 2)).ljust(8) + " units time steps")

  # Plots
  plt.figure("example")
  f, ax = plt.subplots(2, sharex = True)
  ax[0].plot(dates, X)
  ax[0].set_title("Raw signal")
  ax[0].set_ylabel("10$^3$ km$^3$")
  
  ax[1].plot(dates, noise)
  ax[1].set_title("Estimated noise term")
  ax[1].set_ylabel("10$^3$ km$^3$")

  figname = "example.png"
  plt.savefig(figname)
  print("Figure " + figname + " printed")


# If the script is called as is, just run a test-case
if __name__ == '__main__':
  import matplotlib.pyplot as plt
  plt.close("all")
  example()
 
