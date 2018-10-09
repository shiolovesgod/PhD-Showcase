newImg = cell(1,180);

for i = 1:180
    newImg{i}=imrotate(shutDown,i,'crop');
    F(i)= getframe;
end

figure
% while j<3
    tic
    for i =1:60
        imshow(allImages{i})
        drawnow
        F(i)= getframe;
    end
    toc
%     j = input('\n\nWhat is j?\n');
% end

figure
a = 0.5; b=(1-a)/2;
hButton = uicontrol('Style','pushbutton','Units','normalized',...
    'Position',[b,b,a,a]);

   tic
    for i =1:180
        set(hButton,'CData',newImg{i})
        drawnow
    end
    toc

allImages=cell(1,60);
for i =1:60
    
    if i<10
        fname = ['image_f0',num2str(i),'.png'];
    else
        fname = ['image_f',num2str(i),'.png'];
    end
    
        allImages{i} = imread(fname);
end

figure(1)
quitt = false;
hImg = imshow(imread('image_f01.png'));
set(hImg, 'HitTest','off')
set(gca,'ButtonDownFcn',@sampleCallback,'HitTest','off')
set(1,'WindowButtonDownFcn',@(varargin) disp('You Clicked the Figure'))
set(hImg,'WindowButtonDownFcn',@(varargin) disp('You Clicked the Image'))
set(gca,'ButtonDownFcn',@(varargin) disp('You Clicked the Axes'))
k = 0;
while ~quitt && k<3
    for i =1:60
        j = mod(i,60)+1;
        imshow(allImages{j})
        drawnow
        getframe;
        
        if quitt
            break
        end
    end
    k=k+1;
end
