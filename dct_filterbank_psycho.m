

function [y, bitrate]=dct_filterbank_psycho(x,Fs,M,MASK_dB)


bool_plot=0; %set to 1 to plot for each frame

if bool_plot==1
    fprintf('\nPlotting per frame is on.\n');    
else
    fprintf('\nPlotting PER FRAME is OFF! Set bool_plot=1 to see plots.\n');
end

sig=x';
len=length(sig);% signal length
win_length=2^(ceil(log2(20*(Fs/1000))));% frame length
han_win=hanning(win_length)';
sig=[sig zeros(1,win_length-mod(len,win_length))];% pad signal with zeros
len=length(sig);% update the signal length
win_N=ceil(len/win_length);% frame numbers
spectrogram(sig,1024,512,1024,Fs,'yaxis');
title('Input Signal spectrogram');

fb=mdct_filterbank(M);% get the filter banks
inv_fb=fliplr(fb);% get the inversed filter bank
ds_n=M;% downsample step size

y=zeros(M,length(sig));% subband signal
y_ds=zeros(M,ceil(length(sig)/(ds_n)));% downsampled subband signal
y_intp=zeros(M,length(sig));% interpolatation signal    
y_inv=zeros(1,len);
y_q=zeros(size(y_ds));

for i=1:M
    y(i,:)=filter(fb(i,:),1.0,sig);% filter original signal with filterbank
    y_ds(i,:)=y(i,1:ds_n:end);% downsample signal
    y_intp(i,1:ds_n:length(sig))=y_ds(i,:);% interpolation
    y_inv=y_inv+filter(inv_fb(i,:),1.0,y_intp(i,:));% signal reconstruction
end

y_inv=[y_inv(M*2:end), zeros(1,M*2-1)];%adjust for delay


if bool_plot==1
    
figure;
subplot(311);
plot(sig);title('Original signal');

subplot(312);
plot(y_inv);title('Reconstructed signal with error (no quantization)');
hold on;
plot(sig-y_inv);

subplot(313);
plot(sig-y_inv,'r');
title('Error (no quantization)');
end

pause(0.5);


NDFT=win_length;
si=1; ei=NDFT;
win=0;
X_dft=zeros(win_length,win_N);% DFT of windowed frame
SPL=zeros(win_N,win_length/2);% SPL--sound presure level initalization
SPL_band=zeros(M,win_N);%max SPL for each subband
SMR=zeros(size(SPL_band));
maskThr_current=zeros(M,win_N); % threshold of current band
maskThr_prev=zeros(M,win_N); % threshold of previous band
th_mask=zeros(M,win_N); % temporal masking
maxlevel1k=max(abs(fft(han_win.*sin(2*pi*1000*[1:NDFT]/Fs)))); % get maxlevel1k
%maxlevel1k=max(abs(real(fft(sin(2*pi*1000*linspace(0,1/Fs,NDFT)))))); % get maxlevel1k
bits_for_subband=zeros(M,win_N);
bits_for_subband(bits_for_subband<0)=0;
matrix_of_bits=zeros(M,win_N);
samples_in_subband=win_length;


% get the threshold in quiet
f=(Fs/(2000*M))*((2*[[1:M]]-1)/2);
tiq= 3.64*f.^(-0.8) - 6.5*exp(-0.6*(f-3.3).^2) + (10^-3)*f.^4;


while(ei < length(sig)) 
    win=win+1;
    
    %SPL calculation
    time_index_fullband = si:ei;
    X_dft(:,win)=fft(han_win.*sig(time_index_fullband));% DFT of each frame
    SPL(win,:)= 96 + 20*log10(abs(X_dft(1:NDFT/2,win)/maxlevel1k)); % SPL of each frame
    
    
    idx=1:8:NDFT/2;
    dftIndices = zeros(numel(idx),8);                       
        for i = 1:numel(idx)
            dftIndices(i,:)=idx(i):idx(i)+7;
        end
    dftIndices=dftIndices';
    
    
        for bandIndex=1:M
            SPL_band(bandIndex,win)=max(0,max(SPL(win,dftIndices(:,bandIndex))));
        end
    
        
       
    maskThr_current(:,win)=max(conv(SPL_band(:,win),[0.05 0.6 0.3 0.05],'same')-MASK_dB,tiq');
        if win==1
             th_mask(:,win)=max(0, maskThr_current(:,win));
             maskThr_prev(:,win)=zeros(M,1);
        else
            th_mask(:,win)=max(th_mask(:,win-1)*0.9, maskThr_current(:,win));
             maskThr_prev(:,win)=maskThr_current(:,win-1)*0.9;
        end
    
        
      time_index_subband = max(int32(si)/M,1):ei/M; 
    %fprintf('\n%i--%i\n',max(int32(si)/M,1),ei/M);
    
    for ind_subband = 1:M
        
        SMR(ind_subband,win)=SPL_band(ind_subband,win)-th_mask(ind_subband,win);
        bits=SMR(ind_subband,win)/6.02;
        
        if bits>0
        bits_for_subband(ind_subband,win)=round(bits);
        
        else
            bits_for_subband(ind_subband,win)=0;
            
        end
        
        if bits_for_subband(ind_subband,win)>0
           
       y_q(ind_subband,time_index_subband)=myquantizer(y_ds(ind_subband,time_index_subband),bits_for_subband(ind_subband,win));
        end
    
        matrix_of_bits(ind_subband,win)=bits_for_subband(ind_subband,win)*samples_in_subband;
        
    end   
        
    %plots
    if bool_plot==1
        
      
      figure(1)%  subplot(3,1,1);
            plot(1:M,tiq','r');
            hold on;
            plot(1:M,SPL_band(:,win),'b');
            
            plot(1:M,maskThr_current(:,win)','g');
            hold off;
            title(['Frame ',num2str(win)]); 
            legend('Threshold in Quiet','SPL(band)','Masking Threshold');
            xlabel('Band Index');
            ylabel('(in dB)');
        
        figure(2)%subplot(3,1,2);
            plot(1:M,maskThr_current(:,win),'bo');
            hold on;
            plot(1:M,maskThr_prev(:,win),'ro');
            plot(1:M,th_mask(:,win),'k.');
            hold off;
            title(['Frame ',num2str(win)]); 
            legend('Masking Threshold of current window','Masking Threshold of old window','Joint mask');
            xlabel('sub-and index');
            ylabel('Amplitude (dB)');
            pause(0.001)
            
        figure(3) %subplot(3,1,3);
            plotyy(1:length(SMR(:,win)),SMR(:,win),1:length(bits_for_subband(:,win)),bits_for_subband(:,win))
            legend('SMR (dB)','Bits allocated');
            title(['Frame ',num2str(win)]); 
            
        drawnow;

        pause(0.001);
    end
    
    
%next loop indices:
    si=win*win_length+1;
    ei=(win+1)*win_length;
end

sm=sum(matrix_of_bits);
bitrate=sum(sum(matrix_of_bits))/len/M

%Signal reconstruction after quantization

y_inv=zeros(1,len);

for i=1:M
    y_intp(i,1:ds_n:length(sig))=y_q(i,:);% interpolation
    y_inv=y_inv+filter(inv_fb(i,:),1.0,y_intp(i,:));% signal reconstruction
   
end


figure;
y_inv=[y_inv(M*2:end), zeros(1,M*2-1)];


subplot(2,1,1);
plot(sig);
title('Original signal');

subplot(2,1,2);
plot(y_inv)
hold on;
plot(sig-y_inv);
title('Reconstructed with quantization')
soundsc(y_inv,Fs)

figure;    


%for i=1:win_N
    imagesc(linspace(0,len/Fs,512),linspace(0,Fs/2,512),SMR);
    set(gca,'Ydir','Normal');
    colorbar();
    xlabel('Time(s)');
    ylabel('Frequency(Hz)');
    title('SMR');

end
