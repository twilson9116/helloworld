function save_all_figs(save_location)
%%saves all open figures into the specified foler. 
cd(save_location)
h=get(0,'children');
for i=1:length(h);
    saveas(h(i),['figure',num2str(i)],'fig');
end
cd('../')

