function [ cycle ] = f_cycle( adj_N, vect ,k )
%Function to find all the cycles starting and ending at 'k' based on the
%information given in adjacency matrix adj_N,
%vect is the current loop being calculated 
cycle = [];
    for i = 1: size(adj_N,1)
        l_v= size(vect,2);
        if(adj_N(vect(l_v),i)==1 )
            if(size(find(vect==i),2)==0)
                c1 = f_cycle(adj_N,[vect i],k);
            elseif(i==k && vect(l_v-1)~=k)
                c1 = [vect k];
            else
                c1 = [];
            end
            if(size(c1,1)==1 && size(c1,2) ~=0)
                c1 = [c1 zeros(1,size(adj_N,1)+1-size(c1,2))];
            end
            cycle = [cycle;c1];
        end
        
    end


