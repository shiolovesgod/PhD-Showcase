fps = 12;

for percent = 1:99
    barPlot = create();
    barPlot(:,:,2) = 1; barPlot(:,end-percent:end,1) = 1; barPlot(:,end-percent:end,2) = 0; repmat(barPlot, [3,1]);
    imshow(ans);
    pause(1/fps);
end