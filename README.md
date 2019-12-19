# Gaussian process simulator demo in R and Shiny
This is **a simple demonstration of how** (one-dimensional) **Gaussian processes can be simulated**. It is a Shiny webapp with some controls that can draw trajectories of a few Gaussian processes (see the list of kernels below) and calculate a few basic properties of them. This code was shared in hope that someone might benefit from it - eg. as a little program to show off in a class, as a visualization of these processes or just some programming tips to speed up their own way of dealing with these processes.

## Running this program
To run this Shiny app on your device, you'll need to have **R**, **RStudio** and the following packages installed:

* `shiny` (RStudio will ask you to download it upon starting the app for the first time)
* `ggplot2` (used for visualizations)
* `data.table`
* `dplyr`
* `e1071` (used as a faster way to simulate Brownian bridge)

You can install them using RStudio's *Packages* tab or by running the following command: 

```R
install.packages(c('shiny', 'ggplot2', 'data.table', 'dplyr', 'e1071'))
```

Then, to run the app:

1. Download it (either as a zip (the green button to the top left on this page)), or if you're familiar with `git`, you can use `git clone`.
2. (If you've downloaded it as a zip file) Extract it to *any* folder
3. Open the R project in RStudio (either by double clicking the `GaussianProcess-Shiny.Rproj` file or by going to `File -> Open Project...` and pointing to the previously mentioned file)
4. Inside RStudio, open `app.R` (you can use the *Files* tab for this)
5. To the top left corner of RStudio's file editor, click **Run app**

## Supported Gaussian processes
The following **one-dimensional** processes can be drawn:

* Random planes
* Brownian motion [needs axis minimum to be set to 0]
* Brownian bridge [needs axis to be set to [0, 1]]
* Square exponential
* Ornstein-Uhlenbeck [needs axis minimum to be set to 0]
* Periodic

## Metrics calculated by this app
For all the processes

* mean time to reach a certain boundary (like +1 and -1)
* histogram and density of the distribution of reaching these boundaries

For the Brownian bridge

* distribution of maximum of the trajectory's abolute value (how far the process went from 0): histogram and density

## Caveats
Please be aware that this implementation, while it uses a few simple optimizations (mostly: it utilizes matrix multiplications for faster covariance calculations), might still need further work to operate at faster speeds. **It was made for demonstration purposes only and should not be used as-is in a production-ready environment.**

## Resources
If you're interested in Gaussian processes, the following resources may be of help:

* [mathematicalmonk's videos on Gaussian processes](https://www.youtube.com/watch?v=vU6AiEYED9E)
* [the Wikipedia page of Gaussian processes](https://en.wikipedia.org/wiki/Gaussian_process)
