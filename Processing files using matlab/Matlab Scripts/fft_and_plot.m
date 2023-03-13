function [fft_out,q_out]=fft_and_plot(data_line,x_range,plot_flag)

    L=length(data_line);
    pixelsize=x_range/(L-1); %T=1/Fs in the examples
    %x_series=linspace(0,x_range,L);
    x_series=(0:L-1)*pixelsize;
    
    
    y_fft=fft(data_line);
    abs_fft=abs(y_fft);
    x_fft=(1/L)*(1/pixelsize)*(0:L/2);
    %x_fft=(1/(2*pi*range))*(1/pixelsize)*(0:range);
    %f=(1/pixelsize)*(0:L/2)/L;
    

    if exist('plot_flag','var')
        figure;    
        plot(x_fft(2:50),abs_fft(2:50)); 
        xlabel('q (um-1)');
    end
    fft_out=y_fft;
    q_out=x_fft;

 
