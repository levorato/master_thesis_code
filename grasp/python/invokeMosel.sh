for file in *.g; do
        cp -f "$file" input.mos
        mosel -s -c "exec correlationCC.mos"
        cp resultado.dat "${file%.g}-result.dat"
done


