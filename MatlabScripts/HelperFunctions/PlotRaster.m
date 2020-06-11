function PlotRaster(Timestamps,Reference,Window)
%timestamps are the times of the events (spikes usually)
%reference is the event times that the raster is aligned
%window is the range in which events are plotted relative to reference

hold on;
trials=length(Reference);
for trial=1:trials
    times=Timestamps-Reference(trial);
    times=times(times>Window(1) & times<Window(2));
    for stamp=1:length(times)
        plot([times(stamp) times(stamp)],[7+10*(trial-1) 13+10*(trial-1)],'color','k');
    end
end

yticks([100:100:10*trials]);
yticklabels([10:10:trials]);
ylabel('Reference No.');
ylim([0 10*trials+10]);
xlabel('Seconds from reference');
set(gca,'ydir','reverse');