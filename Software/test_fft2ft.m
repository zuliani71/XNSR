% Regression tests for fft2ft frequency axis and one-sided normalization.
Fs_test = 100;

N_even = 16;
OUT_even = fft2ft(fft(randn(N_even,1)),Fs_test);
assert(abs(OUT_even(end,1)-Fs_test/2) < 1e-12);

N_odd = 15;
OUT_odd = fft2ft(fft(randn(N_odd,1)),Fs_test);
f_odd_expected = floor(N_odd/2)*Fs_test/N_odd;
assert(abs(OUT_odd(end,1)-f_odd_expected) < 1e-12);

OUT_matrix = fft2ft(fft(randn(N_even,3)),Fs_test);
f_matrix = OUT_matrix(:,1,1);
assert(abs(f_matrix(end)-Fs_test/2) < 1e-12);
assert(all(diff(f_matrix) > 0));
assert(abs((f_matrix(2)-f_matrix(1))-Fs_test/N_even) < 1e-12);

n_even = (0:N_even-1)';
x_even = 3 + 2*cos(2*pi*2*n_even/N_even) + 4*cos(pi*n_even);
OUT_amp_even = fft2ft(fft(x_even),Fs_test);
assert(abs(abs(OUT_amp_even(1,2))-3) < 1e-12);
assert(abs(abs(OUT_amp_even(3,2))-2) < 1e-12);
assert(abs(abs(OUT_amp_even(end,2))-4) < 1e-12);

n_odd = (0:N_odd-1)';
x_odd = 3 + 2*cos(2*pi*2*n_odd/N_odd);
OUT_amp_odd = fft2ft(fft(x_odd),Fs_test);
assert(abs(abs(OUT_amp_odd(1,2))-3) < 1e-12);
assert(abs(abs(OUT_amp_odd(3,2))-2) < 1e-12);

OUT_amp_matrix = fft2ft(fft([x_even,2*x_even]),Fs_test);
assert(abs(abs(OUT_amp_matrix(1,2,2))-6) < 1e-12);
assert(abs(abs(OUT_amp_matrix(3,2,2))-4) < 1e-12);
assert(abs(abs(OUT_amp_matrix(end,2,2))-8) < 1e-12);

fprintf('N pari: %.12f Hz\n',OUT_even(end,1));
fprintf('N dispari: %.12f Hz\n',OUT_odd(end,1));
fprintf('Ampiezze N pari: DC=%.12f, interna=%.12f, Nyquist=%.12f\n', ...
    abs(OUT_amp_even(1,2)),abs(OUT_amp_even(3,2)), ...
    abs(OUT_amp_even(end,2)));
fprintf('Ampiezze N dispari: DC=%.12f, interna=%.12f\n', ...
    abs(OUT_amp_odd(1,2)),abs(OUT_amp_odd(3,2)));
disp('TUTTI I TEST FFT2FT SONO SUPERATI');
