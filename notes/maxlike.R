Using this information, for any historic or future LTDR (5km) unit of observation, the probability each unit is forest or non-forest is calculated based on a maximum likelihood decision rule.
The probability of each unit containing a forest is calculated, and all units which have a higher probability of inclusion into the "forest" class than "non-forest" class are used in later steps of the analysis.
The estimated probability density function to estimate inclusion into the forested class is calculated using the following equation:
  \begin{equation}
\rho (x | c_{i}) = \frac{1}{(2\pi)^{1/2}*\sigma_{i}} * exp(-\frac{1}{2}*\frac{(x-\mu_{i})^{2}}{\sigma_{i}^{2}})
\end{equation}
where \textit{x} is the LTDR value of the unit being classified as forest or non-forest (represented as classes \textit{c}), \begin{math}\mu_{i}\end{math} is the estimated mean of all LTDR values determined to be included in class \textit{i} and \begin{math}\sigma_{i}^{2}\end{math} is the estimated variance of each class, i.e. \begin{math}x \in forest (c_{1}) \end{math}when:
  \begin{equation}
\rho (x | c_{1}) * \rho(c_{1}) \geq \rho (x | c_{0}) * \rho(c_{0})
\end{equation}
where prior probabilities \begin{math} \rho(c_{i})\end{math} are established based on the probability of a LTDR pixel being classified as forest or non-forest in the year 2000 (i.e., the percent of the landscape occupied by each class).
\par