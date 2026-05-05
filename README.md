
# Final Project Overview

## April 27: Example presentation from 380: https://github.com/EAES-380/final-project-biogeo-biochem.git

Choose a dataset and:  

  * apply data science to explore and wrangle
  * pose 2-3 good questions based on your data exploration
  * answer the questions using methods from EAES 480
  * summarize your data, method, and results in README.md and a 5-min presentation.

## Datasets
  
**Undergrads**  
Daytonna: Help search on 4/8.  
Eli: Help search on 4/8.  
Dylan J.: Has data already.  
Rachel: Help search on 4/8.  
AnaBelle: NA  
Indigo: NA  
Samarth: NA  
Arianna: Discuss with Gavin.  
  
**Grads**  
Aaron: Has data already.  
Amelia: [Forest Inventory Analysis Data for Southeast Alaska](https://github.com/McNicol-Lab/Yellow-cedar-data-compilation)  
Dylan E.: Discuss with Gavin.  
Louis: NA  
Santiago: Discuss with Gavin.  
Wen: [Coastal temperate rainforest soil carbon](https://datadryad.org/dataset/doi:10.5061/dryad.5jf6j1r) [use 2024 version]  

## Good Data Sources

### 1. Our World in Data (OWID) – Climate, Emissions, Energy, and Environmental Indicators  

Link: [www.ourworldindata.org](www.ourworldindata.org) (explore via charts or the CO₂ & Greenhouse Gas Emissions explorer and similar pages for temperature, land use, etc.)  
  
What you get: Clean, country-year (or global) tabular CSVs with long time series (often 1800s–present) on CO₂/GHG emissions (by fuel, per capita, consumption-based), surface temperature anomalies, energy mix, land-use change, biodiversity indicators, and more. Multiple related variables in one dataset for multivariate analysis.

### 2. USGS Earthquake Catalog

Link: [www.earthquake.usgs.gov/earthquakes/search](www.earthquake.usgs.gov/earthquakes/search) (or the real-time feeds at earthquake.usgs.gov/earthquakes/feed/v1.0/csv.php for starters)  
  
What you get: Tabular CSV with columns for time, latitude/longitude, depth, magnitude, location name, type (earthquake vs. explosion), and more. Queryable for global or regional events over any time period (e.g., all M4+ worldwide since 1900).

### 3. Climate TRACE Global Emissions Inventory

Link: [www.climatetrace.org/data](www.climatetrace.org/data) (downloads section)  
  
What you get: ZIP packages of multiple CSVs with country-level or sector-level greenhouse gas emissions (CO₂e by gas, subsector, source) from 2015–2025 (including monthly and projected data). Covers power, manufacturing, fossil fuels, agriculture, waste, etc., with confidence intervals
  
## Deliverables

1.  Proposal - due Wed. April 15, at midnight. (complete proposal.Rmd in your fork)
2.  Presentation - due Wed. May 6, 1-3pm (recording or live via Zoom)
4.  Executive summary - due Sunday, May 10, at midnight. (completed README.md in your fork)


### Proposal

This is a draft of the introduction section of your project as well as a data analysis plan and your dataset.

-   **Section 1 - Introduction:** The introduction should introduce your general

    research question and your data (where it came from, how it was collected,

    what are the cases, what are the variables, etc.).

-   **Section 2 - Data:** Place your data in the \`/data\` folder, and add dimensions and codebook to the README in that folder.
    Then print out the output of and codebook to the README in that folder.
    Then print out the output of `glimpse()` or `skim()` of your data frame.

-   **Section 3 - Data analysis plan:**

    -   The outcome (response, Y) and predictor (explanatory, X) variables you will use to answer your question.

    -   The comparison groups you will use, if applicable.

    -   Very preliminary exploratory data analysis, including some summary statistics and visualizations, along with some explanation on how they help you learn more about your data.
        (You can add to these later as you work on your project.)

    -   The method(s) that you believe will be useful in answering your question(s).
        (You can update these later as you work on your project.)

    -   What results from these specific statistical methods are needed to support your hypothesized answer?

Each section should be no more than 1 page (excluding figures).
You can check a print preview to confirm length.

The grading scheme for the project proposal is as follows.
Note that after you receive feedback for your proposal you can improve it based on the feedback and re-submit it.
If you re-submit, your final score for the proposal will be the average of two scores you receive (first and second submission).

| Total                                | 20 pts |
|--------------------------------------|--------|
| Data description                     | 8 pts  |
| Proposal                             | 10 pts |
| Workflow, organization, code quality | 2 pt   |

### Executive summary

Along with your presentation slides, we want you to provide a brief summary of your project in the README of your repository.

This executive summary should provide information on the dataset you're using, your research question(s), your methodology, and your findings.

The executive summary is worth 30 points and will be evaluated based on whether it follows guidance and whether it's concise but detailed enough.

### Presentation

5 minutes max.
You can either present live during your workshop or pre-record and submit your video to be played during the workshop.

Prepare a slide deck using the template in your repo.
This template uses a package called `xaringan`, and allows you to make presentation slides using R Markdown syntax.
There isn't a limit to how many slides you can use, just a time limit (5 minutes total).
Your presentation should not just be an account of everything you tried ("then we did this, then we did this, etc."), instead it should convey what choices you made, and why, and what you found.

Before you finalize your presentation, make sure your chunks are turned off with `echo = FALSE`.

Presentations will take place during the last workshop of the semester.
You can choose to do your presentation live or pre-record it.
During your workshop you will watch presentations from other teams in your workshop and provide feedback in the form of peer evaluations.
The presentation line-up will be generated randomly.

The grading scheme for the presentation is as follows:

| Total                                                                                                                                                                                                       | 50 pts |
|--------------------------------------------------------|---------|
| Content: Quality of question and match with data       | 10 pts  |
| Methods: Did you use and describe stats accurately?    | 10 pts  |
| Creativity: Did you draw on domain knowledge?          | 10 pts  |
| Time management: Did you use the time well?            | 5 pts   |
| Professionalism: How well did you communciate?         | 5 pts   |
| Organization: Did you tell at a story?                 | 5 pts   |
| Slides: Was use of visuals effective and clear?        | 5 pts   |

### Repo organization

The following folders and files in your project repository:

-   `presentation.Rmd` + `presentation.html`: Your presentation slides
-   `README.Rmd` + `README.md`: Your write-up
-   `/data`: Your dataset in CSV or RDS format and your data dictionary
-   `/proposal`: Your project proposal

Style and format does count for this assignment, so please take the time to make sure everything looks good and your data and code are properly formatted.

## Marking

| Total                            | 100 pts |
|----------------------------------|---------|
| Proposal                         | 20 pts  |
| Presentation                     | 50 pts  |
| Executive summary                | 30 pts  |
