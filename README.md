# Social-and-demographic-predictors-of-survival-in-female-blue-monkeys
Data and code for paper: "Stronger social bonds do not always predict greater longevity in a gregarious primate"
https://onlinelibrary.wiley.com/doi/abs/10.1002/ece3.3781

Survival analysis on 83 adult female blue monkeys, assessing how social factors influence risk of death.
Analysis of both fixed-time (multi-year) and time-varying (annual) social correlates of survival.

Variables in data set include:
**strength.3 & strength.6**: Average Dyadic Sociality Index (DSI based on grooming and spatial association) with top 3 and top 6 partners\
**cons.3 & cons.6** : % of top partners consistent from year to year\
**st.co3 & st.co6** : Strength consistency class. Levels:
  1 = above average strength & above average consistency
  2 = above average strength & below average consistency
  3 = below average strength & above average consistency
  4 = below average strength & below average consistency
 **rank** : Standardized dominance rank - test effect of contest competition
 **af.groupmates**: Average number of adult female groupmates - test effect of scramble competition
 **age.first.rep** : Age at first reproduction in days - test effect life history trade off
 **period.pres** : Number of years present
 **date.of.birth.estimate.range.years** : Range of error in birth date estimate
 
 Analyses in script include:
 1. Repeatability of annual social and demographic variables
 2. Fixed-time effect Cox models
    Testing quadratic term of number of female group mates
 3. Time-varying covariate Cox models
 4. Compare estimated coefficients to distribution of coefficients based on 1000 random permtuations of underlying social matrices
