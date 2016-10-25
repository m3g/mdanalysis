!
! Program blockavg: Computes the estimated error of a correlated
!                   time series of data using the "Blocking Average" method
!                   described in:
!
!                   Flyvbjerg, Petersen, Error estimates on averages of correlated data.
!                   J. Chem. Phys. 91, 461 (1989)
!

program blockavg


  






end program blockavg


Estou enviando este e-mail, originalmente para a Emília,
para todos, porque acho que é o do interesse de todo mundo.
Trata-se do problema de estimar qual é o erro na medida
da média de uma propriedade que flutua, quando temos uma ou mais
simulações.

Oi Emília,

Estive pensando um pouco nos erros das médias, e a médida
que a gente fez ainda não é exatamente correta.

Cada simulação tem sua média, e as flutuações dessa média
dependem do tempo de cada simulação. Por exemplo, se cada
uma das simulações fosse muito longa, todas deveriam dar
a mesma média ou, dito com um pouco mais de cuidado, as
flutuações da média seriam muito pequenas.

Portanto, as flutuações da média dependem do tamanho das
simulações. Isto quer dizer que, dado um tamanho de simulação
você tem um perfil característico das "flutuações" da sua média.
Ou seja, se suas simulações são de 10 ns, você vai ter
distribuição de médias com uma certa dispersão, enquanto que
se suas simulações forem de 100 ns, você vai ter uma distribuição
de médias com uma dispersão menor.

No entanto, se fizermos 10 vezes mais simulações de 10 ns
em relação a simulações de 100 ns, deveríamos ser capazes
de prever a média com a mesma precisão, já que o tempo total
de amostragem é o mesmo.

Ou seja, se suas simulações são de 10 ns, mas você faz mais
simulações, a qualidade da sua estimativa da média deveria
melhorar. No entanto, a distribuição das médias é independente
do número de simulações de 10 ns que você faz, ela é uma
característica desse tempo de simulação. No caso extremo,
em que cada simulação tem um único passo, a distribuição
das "médias" é igual à distribuição dos valores possíveis
da variável de interesse. Como sabemos, a distribuição
dos valores não é um erro, é uma reflexo da variabilidade
intrínseca dos valores no sistema, portanto não podemos tomar
o desvio padrão das médias como o erro na estimativa da
média.

Bom, esse problema já foi estudado e já foi resolvido. O problema
se chama "erro da estimativa da média". O erro na estimativa
da média, quando cada medida é totalmente descorrelacionada
da outra medida, é (usando a variância como medida de erro):

variância_da_estimativa_da_média =
     variância_da_distribuição_da_propriedade / sqrt( número_de_médias )

( desvio padrão é a raiz quadrada da variância).

Ou seja, seu erro é menor quanto maior o número de médias que
você faz. No seu caso, o procedimento é o seguinte: Calcular
o desvio padrão da propriedade (número de h-bonds), usando
todas as simulações juntas. Esse é o desvio padrão da distribuição
dos valores possíveis do número de h-bonds.

Se você fez 10 simulações, e considerando que cada simulação
está descorrelacionada da outra, o erro na estimativa da média
total vai ser dado pela conta acima, ou seja, tem que dividir
pela raiz quadrada do número de simulações.

A confiança no resultado é de 68% para mais ou menos
um desvio padrão e 95% para mais ou menos dois desvios padrão.
(existe um teorema, chamado "teorema do limite central", que
diz que a distribuição das estimativas das médias é gaussiana
independentemente da natureza da distribuição da propriedade,
por isso podemos usar as propriedades das distribuições normais
sem medo).

Essa conta é simples no caso em que assumimos que
as simulações independentes estão descorrelacionadas,
ou seja, que o tempo de correlação da variável medida é muito
menor que o tempo de cada simulação individual. Se não for
assim, ou se tivermos uma única simulação (por exemplo,
dez vezes mais longa), o método bom para fazer a mesma
coisa se chama "block averaging". Está descrito em vários
lugares, por exemplo aqui:

http://sachinashanbhag.blogspot.com.br/2013/08/block-averaging-estimating-uncertainty.html

(o artigo importante é Flyvbjerg e Petersen, J. Chem. Phys. 91, 461 (1989).
