pkg load control
pkg load statistics

p = 18/37
N = 2*10**3
xGrid = linspace(1,7,7)

function x = rouletteSpins(p)
  x = 0
  while true
    x = x + 1
    if rand() < p
      return
    endif
  endwhile
endfunction
mc = []
for i = 1:N
  mc(end+1) = rouletteSpins(p)
endfor

mcEstimate = []
for i = xGrid
  mcEstimate(end+1) = sum(mc==i)/N
endfor

gPmf = geopdf(xGrid-1,p)

stem(xGrid, mcEstimate)
hold on
stem(xGrid,gPmf,"x")
xlabel("X")
ylabel("Probability")
legend("MC estimate","PMF")
hold off