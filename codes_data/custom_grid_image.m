function custom_grid_image(xgrid,ygrid,varargin)

if nargin==2
    va=[];
else
    va=varargin;
end
start_hold=ishold;
hold on
if ~isempty(xgrid)
    for i=1:length(xgrid)
        h=line('xdata',xgrid(i)*[1 1],'ydata',ylim);
        if ~isempty(va)
            set(h,va{:})
        end
    end
end
if ~isempty(ygrid)
    for i=1:length(ygrid)
        h=line('xdata',xlim,'ydata',ygrid(i)*[1 1]);
        if ~isempty(va)
            set(h,va{:})
        end
    end
end
if ~start_hold
    hold off
end

end