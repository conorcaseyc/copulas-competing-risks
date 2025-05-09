# Using Copulas to Model Prepayment and Default in Competing Risks Analysis

This repository contains the code accompanying the paper **"Using Copulas to Model Prepayment and Default in Competing Risks Analysis"**, produced as part of the STU44003 Data Analytics module at Trinity College Dublin.

## Overview
This project explores the application of copula models to analyze the dependency between competing financial risks, specifically loan prepayment and default. The approach combines flexible marginal modelling (e.g., Weibull, gamma, and mixture distributions) with Archimedean copulas (Frank, Clayton, Gumbel) to construct joint distributions of event times. Both simulation studies and real-world analysis on Lending Club data are included.

### Highlights
- Simulated competing risks data using **Frank** and **Clayton** copulas
- Fitted marginal distributions (Gamma, Exponential, Weibull, and Mixture Models)
- Joint likelihood construction and **maximum likelihood estimation (MLE)**
- **Bootstrap procedures** to assess parameter uncertainty
- Real-world application to over 2 million Lending Club loans
- Comparative evaluation of **Gumbel**, **Frank**, and **Clayton** copulas