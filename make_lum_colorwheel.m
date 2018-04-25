

hsl_colorwheel = NaN(size(newest_colorwheel));

for ii = 1:length(newest_colorwheel)

    hsl_colorwheel(ii,:) = rgb2hsl(newest_colorwheel(ii,:)./255);
    
    hsl_colorwheel(ii,1) = hsl_colorwheel(ii,1)*360;
    hsl_colorwheel(ii,2:3) = hsl_colorwheel(ii,2:3).*100;

end

hsl_colorwheel(1:360,3) = 50;



lum_colorwheel = NaN(size(hsl_colorwheel));

for ii = 1:length(hsl_colorwheel)

    lum_colorwheel(ii,1) = hsl_colorwheel(ii,1)/360;
    lum_colorwheel(ii,2:3) = hsl_colorwheel(ii,2:3)./100;
    
    lum_colorwheel(ii,:) = hsl2rgb(lum_colorwheel(ii,:)).*255;

end



