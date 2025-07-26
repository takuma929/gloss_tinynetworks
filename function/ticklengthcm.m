function ticklengthcm(axh,cm)
% Standardize tick mark length across figures (from MATLAB FileExchange and modified)
    p = inputParser;
    p.addRequired('axh', @(x) (isscalar(x) && all(isgraphics(x,'axes') | isgraphics(x,'colorbar'))) || isa(x,'matlab.graphics.Graphics'));
    p.addRequired('cm', @(x) isscalar(x) && x >= 0);
    p.parse(axh, cm);
    assert(~verLessThan('matlab','8.4.0'), 'ticklength only works with MATLAB R2014b or later.');
    for i = 1:numel(axh)
        if isgraphics(axh(i),'axes') || isgraphics(axh(i),'colorbar')
            oriignalunits = axh(i).Units;
            axh(i).Units = 'centimeters';
            pos = axh(i).Position;
            longest = max(pos(3:4));
            newlength = cm / longest;
            if isgraphics(axh(i),'axes')
                axh(i).TickLength = [newlength, axh(i).TickLength(2)];
            elseif isgraphics(axh(i),'colorbar')
                axh(i).TickLength = newlength;
            end
            axh(i).Units = oriignalunits;
        end
    end
end
