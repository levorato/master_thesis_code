#Isso e um comentario, INICIO DO SCRIPT
#o script precisa entrar na pasta onde o programa esta
export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
cd /home_nfs/mlevorato/mestrado/Build_Oscar
#agora o ambiente do mpi precisa ser iniciado. O arquivo de maquinas
#e criado dinamicamente pelo mpi, e seu caminho fica guardado na variavel
#chamada $PBS_NODEFILE.
#para verificar no arquivo de saida como o arquivo de maquinas ficou, podemos
#dar cat no arquivo de maquina
cat $PBS_NODEFILE
# O numero de maquinas que o PBS forneceu pode ser obtido pela contagem do
#numero de linhas (sem repeticao) do arquivo $PBS_NODEFILE
NUMERO_DE_MAQUINAS=`sort $PBS_NODEFILE | uniq | wc -l`
#para verificar que o numero de maquinas esta correto, em relacao ao arquivo de maquinas
cat $NUMERO_DE_MAQUINAS
#agora, iniciar o ambiente do mpi, com o comando mpdboot
mpdboot -n $NUMERO_DE_MAQUINAS -f $PBS_NODEFILE
#e possivel determinar o numero de processos, pelo numero de entradas no arquivo de maquinas
NUMERO_DE_PROCESSOS=`wc -l $PBS_NODEFILE`
#verificar no arquivo de saida, se o numero de processos era o esperado
cat $NUMERO_DE_PROCESSOS
#executar o programa mpi
mpiexec -n $NUMERO_DE_PROCESSOS ./MestradoMario -l 2 --iter=800 --profile=true --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Gisele/Random"
#por ultimo, finalizar o ambiente mpi
mpdallexit
#FIM DO SCRIPT

