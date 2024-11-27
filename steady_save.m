
function steady_save
cur='dip';
    for n=1:4
    name_fig=strcat(cur,num2str(n),".png");
    saveas(figure(n),name_fig)
    end
end