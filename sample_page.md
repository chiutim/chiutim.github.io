## Monte Carlo Simulation with R

<img src="images/rainwater.png?raw=true"/>

### Parameters
**Monthly Rainfall Amounts:**
- Used R package fitdistplus to fit various distributions (normal, Weibull, and exponential) to each month’s yearly totals.
- Normal distribution produced lowest p-value on the Kolmogorov-Smirnov Goodness of Fit test (0.20-0.45), still not significant
- Allow user parameter for how to run monthly rainfall random variable:
  - Normal Distribution with unique mean and SD for each month
  - Weibull Distribution with unique shape and scale for each month
  - Exponential Distribution with unique rate for each month
  - Custom PDF for each month using R’s pdqr library to generate unique distribution to pull random values.
_Note: For normal and custom distributions, negative rainfall amounts are discarded, and a new random amount is generated._

**Simulation Parameters:**
- Number of iterations: 1,000 or user specified
- Number of years per iterations: 30 or user specified
- Roof Area (in sqft): 3,000 or user specified
- Main Tank Volume: 25,000 or user specified
- Minimum/Maximum Monthly Water Use By Rancher: [4,000-5,200]* or user specified
- Starting Tank Level Per Iteration: 10,000 or user specified

**Hard Parameters:**
- Minimum/Maximum Rainfall Capture Efficiency: [0.90-0.98]*

_* Simulation uses uniform continuous distribution between min/max values_


### Algorithm Pseudocode
1. Read in historic monthly rainfall amounts
2. Get user inputs on simulation parameters
3. Fit user-specified distribution to rainfall data for each month (12 distribution models)
4. Generate tracking variables to store simulation data
5. For each iteration
6. Generate overflow and empty counters for current iteration
  - For each year
    - For each month
      - Pull random rain amount from current month’s distribution (if negative, pull again)
      - Calculate total rainfall = roof area * (rain amount / 12) * 7.48052*
      - Remove water captured given random capture rate
      - Calculate random monthly usage by rancher
      - Calculate total difference = monthly rain – monthly usage
      - Add total difference to current tank amount
      - If current tank amount > capacity
        - Remove excess and add one to overflow counter
      - Else if current tank <= 0
        - Set current amount to 0 and add one to empty counter
      - Add all calculations to simulation data tracker
  - Add iteration number of empties and number of overflows to tracking variable


### [R program](/program/montecarlo.r)
