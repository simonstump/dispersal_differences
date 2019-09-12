function out = make_neighbor(dist,SIZE,type)

%out = make_neighbor(dist,SIZE,type)
%
%This function is used for generating a neighborhood matrix.  Basically, if you have a 
%SIZE x SIZE matrix, then the first column will be a list of all of the element a distance 
%"dist" from the first entry, the second column will be a list of the elements a distance %"dist" from the second entry, etc.  The final output will be a (4*dist)xSIZE^2 matrix.
%
%Here, "distance" is defined using Manhattan distances.  Also, if type=1, then we assume 
%there are warping boundary conditions.  If type=2, then we do not (and if something is 
%too close to the edge, it will have fewer neighbors; those "missing" neighbors will be 
%marked with the number SIZE^2+1.
%
%This code was adapted from code used for Stump et al. (2018, Journal of the Royal Society
%interface).


bob=[1:SIZE^2];
bob2=reshape(bob,SIZE,SIZE)';

%Offset defines how far away the neighbors are.  For for example, the neighbor 2 above
%  will be the number -2, and the one two to the right will be +(2*dist).  

%First, I start with dist above, dist below, dist left, and dist right
offset=[0 -dist; 0 dist; dist 0; -dist 0];

for i =1:dist-1

    %Then, I add more.  I add one that is i left and (dist-i) below, one that is i
    %  right and (dist-i) below, i left and (dist-1) above, and i right and (dist-i)
    %  above.

    offset=[offset; i, (-dist+i); -i, (dist-i); -i, (-dist+i); i, (dist-i)];
end


out=zeros(dist*4,SIZE^2);


theX=repmat(1:SIZE,SIZE,1);
theX=reshape(theX,SIZE^2,1);

theY=repmat(1:SIZE,1,SIZE)';

all=[theX theY];

if(type==1)  %if warping boundary conditions
    
    for i=1:SIZE^2
        x1=all(i,1);
        y1=all(i,2);
                
        x2=mod(offset(:,1)+x1,SIZE);
        y2=mod(offset(:,2)+y1,SIZE);
        
        x2(x2==0)=SIZE;
        y2(y2==0)=SIZE;
        
        for j=1:length(x2)
            out(j,i)=bob2(x2(j),y2(j));
        end
        
    end
    
elseif(type==2)  %if not warping boundary conditions
    for i=1:SIZE^2
        x1=all(i,1);
        y1=all(i,2);
        
        xT=offset(:,1)+x1;
        yT=offset(:,2)+y1;
        
        x2=mod(offset(:,1)+x1,SIZE);
        y2=mod(offset(:,2)+y1,SIZE);
        
        x2(x2==0)=SIZE;
        y2(y2==0)=SIZE;
        
        stayIn=(xT==x2).*(yT==y2);
        for j=1:length(x2)
            if(stayIn(j))
                out(j,i)=bob2(x2(j),y2(j));
            else
                out(j,i)=SIZE^2+1;
            end
        end
        

    end
else
    'error, wrong type'
    return
end

%out

