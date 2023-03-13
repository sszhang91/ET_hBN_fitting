function figformatimg(fig_in)
    
    if nargin>0
        set(fig_in,'EdgeColor','none');    
    end
    
    set(gca,'FontSize',13,'LabelFontSizeMultiplier',1,...
        'TitleFontSizeMultiplier',1,'XColor','k','Box','off','TickDir','in');
    axis tight;

end
