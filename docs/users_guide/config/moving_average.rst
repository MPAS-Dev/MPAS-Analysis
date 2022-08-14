.. _config_moving_average:

Moving Average
==============

By default, some time series have a 12-month moving average while others do
not include a moving average (``movingAverageMonths = 1``).  To perform
a moving average (e.g. over 12 months), set::

  movingAverageMonths = 12

This can be useful for taking out the seasonal cycle to better examine annual
mean tends.

