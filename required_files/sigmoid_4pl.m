function [y] = sigmoid_4pl(x, max_asymptote, slope, shift)


  %y =  p(4)+ ((p(1)-p(4))./ (1 + ((x/p(3))*p(2))));
  
  y = max_asymptote./(1+exp(-1*slope.*(x-shift)));
 
  
 
  
    

end

