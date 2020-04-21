function [y,fy]=conn_filter(rt, filter, x, option)

if nargin < 4, option='full'; end

fy=fft(x,[],1);
f=(0:size(x,1)-1);
f=min(f,size(x,1)-f);

switch(lower(option))
    case 'full'
        
        idx=find(f<filter(1)*(rt*size(x,1))|f>=filter(2)*(rt*size(x,1)));
        %idx=idx(idx>1);
        fy(idx,:)=0;
        y=real(ifft(fy,[],1))*2*size(x,1)*(min(.5,filter(2)*rt)-max(0,filter(1)*rt))/max(1,size(x,1)-numel(idx));
        
    case 'partial'
        
        idx=find(f>=filter(1)*(rt*size(x,1))&f<filter(2)*(rt*size(x,1)));
        %if ~any(idx==1), idx=[1,idx]; end
        y=fy(idx,:);
        
end

end