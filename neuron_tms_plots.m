L=1.7061e-05
omega=1.7799e+04
alpha=185.0018
th=[0 10 20 30 40 50 60 70 72 74 76 78 80 81 82 83 84 85 86 87 88 89 90];
thguess=3000/L*[   6.5994    6.6757    6.9809    7.5531    8.5449   10.1471 ...
    12.8937   18.6157   20.5994   22.8882   25.9399 ...
   30.0598   35.7056   39.3677   43.6401   49.4385   56.7627   66.5283  ...
   81.1768  103.1494  141.6016  227.0508 571.2891]
thcorr=1.0e+10 *[  
    0.1193    0.1207    0.1262    0.1365    0.1545    0.1834    0.2331    0.3365    0.3724    0.4138    0.4690 ...
    0.5434    0.6455    0.7117    0.7889    0.8265    0.9490    1.1123    1.3572    1.5842    1.9822    2.2520 ...
    2.3642];
  %%

semilogy(th,thcorr(1:end)/10^6,'LineStyle','-','Marker','none','Color',"#000080",'LineWidth', 3)% blue
hold on
semilogy(th,thguess(1:end)/10^6,'LineStyle','-.','Marker','none','Color',"#FFC000",'LineWidth', 3) %Gold
semilogy(th(end),thcorr(end)/10^6,'LineStyle','none','Marker','o','Color',"#000080",'LineWidth', 3)% blue
xlabel("Angle (degrees)", 'FontSize', 12);
ylabel("Activation Threshold", 'FontSize', 12);
xlim([0 90]);
legend('Bidomain Integral Equation','Cable Equation')