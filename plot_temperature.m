

time_labels = [];

for i=1:20:size(T_plot,2)
    %plot(xc(floor(N/2)-10:floor(N/2))-xc(floor(N/2)), T_plot(floor(N/2)-10:floor(N/2),i) - 273);
    plot(xc(end-20:end-1), T_plot(end-20:end-1,i) - 273);

    hold on
    label = sprintf("Time %d", floor(pplot_times(i)));
    time_labels = [time_labels label];
end

legend(time_labels, 'Location', 'best')

xlabel('Distance from opening (m)')
ylabel('Temperature (C)')

hold off