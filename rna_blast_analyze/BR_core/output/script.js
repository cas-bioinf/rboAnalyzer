<script type="text/javascript">

    var BtnSelectSeqsTxt = "Select all Seqs."
    var BtnUnselectSeqsTxt = "Unselect all Seqs."
    var BtnSelectStrsTxt = "Select all Structs."
    var BtnUnselectStrsTxt = "Unselect all Structs."
    var BtnSortEvalAscTxt = "Sort Eval asc."
    var BtnSortEvalDescTxt = "Sort Eval desc."

    function select_class_from_colect(colect, cl) {
        console.log("function select from collect start");
        var i;
        for (i = 0; i < colect.length; i++) {
            if (colect[i].getAttribute("class") === cl) {
                console.log(colect[i]);
                return colect[i];
            };
        };
    }

    function mypicload(ete) {
        console.log("function mypicload start")
        ete.href = select_class_from_colect(ete.children, "rnapic").getAttribute("src")
    }

    function exportSelected(classname) {
        console.log("export selection")
        /*
        Find all checkboxes
        return data to fasta
        */
        var checkboxes;
        var checked = new Array(0);
        var i;

        checkboxes = document.getElementsByClassName(classname)

        for (i = 0; i < checkboxes.length; i++) {
            if (checkboxes[i].checked) {
                checked.push(checkboxes[i])
            };
        };
        return checked
    }

    function buildFastaStructure(seqname, sequence, structure) {
        return ">" + seqname + "\n" + sequence + "\n" + structure + "\n"
    }

    function buildFasta(seqname, sequence) {
        return ">" + seqname + "\n" + sequence + "\n"
    }

    function writeSelectionStructures() {
        console.log("write selection fasta")
        var sel = exportSelected("individualStructureCheckbox");
        var i;
        var exportdata = "";

        for (i = 0; i < sel.length; i++) {
            var pnodedata = sel[i].parentNode.parentNode.parentNode.parentNode.dataset;
            var structure = sel[i].parentNode.parentNode.dataset.brna_secondary_structure;

            exportdata += buildFastaStructure(pnodedata.brna_seqname, pnodedata.brna_sequence, structure)
        }
        window.open("data:application/txt," + encodeURIComponent(exportdata), "_self");
    }

    function placeSelection(boxlist, val2place) {
        for (i = 0; i < boxlist.length; i++) {
            boxlist[i].checked = val2place;
        };
    }

    function selectAllseqs(button) {
        var checkboxes;
        var i;
        checkboxes = document.getElementsByClassName("individualSequenceCheckbox")
        if (button.textContent === BtnSelectSeqsTxt) {
            button.textContent = BtnUnselectSeqsTxt;
            placeSelection(checkboxes, true);
        } else {
            button.textContent = BtnSelectSeqsTxt;
            placeSelection(checkboxes, false);
        };
    }

    function selectAllstrs(button) {
        var checkboxes;
        var i;
        checkboxes = document.getElementsByClassName("individualStructureCheckbox")

        if (button.textContent === BtnSelectStrsTxt) {
            button.textContent = BtnUnselectStrsTxt;
            placeSelection(checkboxes, true);
        } else {
            button.textContent = BtnSelectStrsTxt;
            placeSelection(checkboxes, false);
        };
    }

    function writeSelectionFasta() {
        console.log("write selection fasta")
        var sel = exportSelected("individualSequenceCheckbox")
        var i;
        var exportdata = "";

        for (i = 0; i < sel.length; i++) {
            var pnodedata = sel[i].parentNode.parentNode.parentNode.dataset;
            exportdata += buildFasta(pnodedata.brna_seqname, pnodedata.brna_sequence)
        }
        window.open("data:application/txt," + encodeURIComponent(exportdata), "_self");
    }

    function GetHits() {
        var divs;
        divs = document.getElementsByClassName("onehit");
        return divs
    }

    function sortByEval(button) {
        var parent = document.querySelector('.hits')

        var a = Array.from(parent.children)

        if (button.textContent === BtnSortEvalDescTxt) {
            a.sort(
                eval_desc
            ).forEach(
                function(ele) {
                    parent.appendChild(ele);
                }
            )
            button.textContent = BtnSortEvalAscTxt;
        } else {
            a.sort(
                eval_asc
            ).forEach(
                function(ele) {
                    parent.appendChild(ele);
                }
            )
            button.textContent = BtnSortEvalDescTxt;
        };
    }

    function eval_desc(a, b) {
        return Number.parseFloat(a.dataset.eval) < Number.parseFloat(b.dataset.eval);
    }

    function eval_asc(a, b) {
        return Number.parseFloat(a.dataset.eval) > Number.parseFloat(b.dataset.eval);
    }
    console.log("script loaded correctly");
</script>
<script type="text/javascript" src="https://www.ncbi.nlm.nih.gov/projects/sviewer/js/sviewer.js"></script>
