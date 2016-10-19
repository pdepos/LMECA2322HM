function [ height_x ] = height(x)
   if x < 0.04560
    end
    height_x = (-1)*sqrt((0.2543^2)-((x-0.03)^2))+0.2573;
    if (0.04554 < x) && (x < 0.1026)
    height_x = 0.0560*(x+0.015); 
    end
    if (0.1026 < x)
    height_x = (1)*sqrt((0.1537^2)-((x-0.12)^2))-0.1462; 
    end
 
end

