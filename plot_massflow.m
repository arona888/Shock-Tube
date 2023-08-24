plot(plot_times(2:end),-diff(totmass)./diff(plot_times))
xlabel('Time (s)')
ylabel('Mass flow (kg/s)')