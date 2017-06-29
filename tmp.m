fig_array=get(0,'children');

mkdir('figures');
cd('figures');
for i=1:length(fig_array);
savefig(fig_array(i),[fig_array(i).Name '.fig'],'compact')
end
 cd('../')
  close all