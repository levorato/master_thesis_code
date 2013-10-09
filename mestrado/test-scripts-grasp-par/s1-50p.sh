#Isso e um comentario, INICIO DO SCRIPT
#o script precisa entrar na pasta onde o programa esta
export LD_LIBRARY_PATH=$HOME/lib:/opt/intel/mpi/3.2.2.006/lib64:$LD_LIBRARY_PATH
cd /home_nfs/mlevorato/graspcc-1.0/build
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
echo $NUMERO_DE_MAQUINAS
#agora, iniciar o ambiente do mpi, com o comando mpdboot
mpdboot -n $NUMERO_DE_MAQUINAS -f $PBS_NODEFILE
#e possivel determinar o numero de processos, pelo numero de entradas no arquivo de maquinas
NUMERO_DE_PROCESSOS=`wc -l $PBS_NODEFILE | cut -f 1 -d' '`
#verificar no arquivo de saida, se o numero de processos era o esperado
echo $NUMERO_DE_PROCESSOS
#executar o programa mpi
mpiexec -n $NUMERO_DE_PROCESSOS ./graspcc -l 1 --iter=50 --alpha=0.8 --input-file-dir "/home_nfs/mlevorato/mestrado/tests/Instances/Social Media/slashdot-undirected/part0" --output-folder "./output/sdot-l1-iter50-0.8-$NUMERO_DE_PROCESSOS-2h" --gain-function-type=0 --time-limit=7200
#por ultimo, finalizar o ambiente mpi
mpdallexit
#FIM DO SCRIPT

