CEA_RESULTS = CEA('prob','rocket','equilibrium', ...    %problem type
    'o/f',4,'p,bar',(3.4e6)./(1e5), ...   % o/f and chamber pressures
    'pi/p', 34, ...                         %chamber/exit pressure ratios
    'reac','fu','RP-1','wt%',100.,'t(k)',298.15, ...    %RP-1 data
    'ox','N2O','wt%',100,'t(k)',300, ...                %N2O data
    'output','transport','end');