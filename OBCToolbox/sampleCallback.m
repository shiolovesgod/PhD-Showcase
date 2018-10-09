function sampleCallback(hObj, eData)
type = get(hObj,'SelectionType');

switch type
    
    case 'normal' %single click 
        disp('You clicked me')
        set(findobj('Tag','animation'),'UserData',true)
    case 'open' %double click
        allImages=cell(1,60);
        for i =1:60
            
            if i<10
                fname = ['image_f0',num2str(i),'.png'];
            else
                fname = ['image_f',num2str(i),'.png'];
            end
            
            allImages{i} = imread(fname);
        end

        quitt = false;
        k = 0;
        while ~quitt && k<3
            for i =1:60
                j = mod(i,60)+1;
                hImg = imshow(allImages{j});
                set(hImg,'Tag','animation','HitTest','off');
                drawnow               
                if get(hImg,'UserData')
                    quitt = true;
                    break
                end
            end
            k=k+1;
        end
end
