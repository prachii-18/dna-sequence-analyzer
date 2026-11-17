#importing all necessary libraries
import streamlit as st
import matplotlib.pyplot as plt

#Browser tab title, icon and centered page layout
st.set_page_config(page_title="DNA Sequence Analyzer", page_icon="🧬", layout="centered")

#Helper Functions
def clean_sequence(raw):
    return "".join(raw.split()).upper() #Remove spaces/newlines and convert to uppercase


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return "".join(complement[base] for base in reversed(seq))


def find_motif(seq, motif):
    positions = []
    start = seq.find(motif)
    while start != -1:
        positions.append(start + 1) #1-based position, overlapping matches are included
        start = seq.find(motif, start + 1)
    return positions


def set_example():
    st.session_state.seq_input = "ATGCGTACGTTAGCCGATATCGCGTAGCTAGGCTA"


# Page Header
st.title("🧬 DNA Sequence Analyzer 🧬")
st.caption("A simple Bioinformatics tool for DNA Sequence Analysis and Visualisation")

# Input Sequence
st.header("Input DNA Sequence")
st.write("Upload a FASTA file or type/paste the sequence below.")

uploaded = st.file_uploader("Upload FASTA file", type=["fasta", "fa"])

if uploaded:
    raw_text = uploaded.read().decode("utf-8-sig", errors="replace")
else:
    st.text_area("Or enter DNA sequence:", key="seq_input", height=120)
    st.button("Use example sequence", on_click=set_example)
    raw_text = st.session_state.get("seq_input", "")

#Strip each line once so headers with leading/trailing spaces are still detected
lines = [line.strip() for line in raw_text.splitlines()]
if sum(line.startswith(">") for line in lines) > 1:
    st.warning("Multiple records found. They are combined into one sequence, so please use one record at a time.")
sequence = clean_sequence("".join(line for line in lines if not line.startswith(">")))

# Analysis
if sequence:
    # Input validation
    invalid = set(sequence) - set("ATGC")
    if invalid:
        st.error(f"Ambiguous or invalid characters found: {', '.join(sorted(invalid))}. "
                 "Please use only A, T, G and C.")
    else:
        # Counts
        counts = {base: sequence.count(base) for base in "ATGC"}
        total_length = len(sequence)
        at = (counts['A'] + counts['T']) / total_length * 100
        gc = (counts['G'] + counts['C']) / total_length * 100

        # Derived values, computed once and reused in the app and the download
        ratio_text = f"{gc / at:.2f}" if at > 0 else "Cannot compute as AT content is zero"
        if at > gc:
            composition = "AT rich"
        elif gc > at:
            composition = "GC rich"
        else:
            composition = "Balanced"
        rev_comp = reverse_complement(sequence)
        is_palindrome = sequence == rev_comp
        rna = sequence.replace("T", "U")

        # Summary
        st.divider()
        st.header("Sequence Summary")
        c1, c2, c3 = st.columns(3)
        c1.metric("Total Bases", total_length)
        c2.metric("GC Content", f"{gc:.2f}%")
        c3.metric("AT Content", f"{at:.2f}%")

        rows = [
            {"Base": b, "Count": c, "Percentage": f"{c / total_length * 100:.2f}%"}
            for b, c in counts.items()
        ]
        st.dataframe(rows, hide_index=True, use_container_width=True)

        # Composition analysis
        st.subheader("Composition Analysis")
        st.write(f"GC/AT Ratio: {ratio_text}")

        if composition == "AT rich":
            st.info("Sequence is AT rich.")
        elif composition == "GC rich":
            st.info("Sequence is GC rich.")
        else:
            st.info("Sequence has balanced AT and GC content.")

        # Reverse complement
        st.subheader("Reverse Complement")
        st.text_area("Reverse Complement", value=rev_comp, height=100)

        if is_palindrome:
            st.info("This sequence is a reverse-complement palindrome")
        else:
            st.info("This sequence is not a reverse-complement palindrome")

        # Transcription
        st.subheader("RNA Transcript")
        st.text_area("RNA", value=rna, height=100)

        # Motif search
        st.subheader("Motif Search")
        motif = clean_sequence(st.text_input("Enter a motif to search for (e.g. GAATTC):"))
        motif_result = ""
        if motif:
            if set(motif) - set("ATGC"):
                st.error("Motif can only contain A, T, G and C.")
            elif len(motif) > total_length:
                st.warning("Motif is longer than the sequence.")
            else:
                positions = find_motif(sequence, motif)
                if positions:
                    st.success(f"Found {len(positions)} match(es) of {motif}.")
                    st.write("Start positions (1-based): " + ", ".join(map(str, positions)))
                    motif_result = f"Motif {motif} found at positions: {', '.join(map(str, positions))}\n"
                else:
                    st.info(f"Motif {motif} was not found in the sequence.")
                    motif_result = f"Motif {motif} not found\n"

        # Download
        base_lines = "\n".join(
            f"{b}: {counts[b]} ({counts[b] / total_length * 100:.2f}%)" for b in "ATGC"
        )
        result_text = (
            f"Length: {total_length}\n"
            f"{base_lines}\n"
            f"GC Content: {gc:.2f}%\n"
            f"AT Content: {at:.2f}%\n"
            f"GC/AT Ratio: {ratio_text}\n"
            f"Composition: {composition}\n"
            f"Reverse-complement palindrome: {'Yes' if is_palindrome else 'No'}\n"
            f"Reverse Complement: {rev_comp}\n"
            f"RNA Transcript: {rna}\n"
            f"{motif_result}"
        )
        st.download_button("Download results", result_text, file_name="dna_analysis.txt")

        # Visualizations
        st.divider()
        st.header("Nucleotide Composition Visualizations")
        colors = ['coral', 'lightseagreen', 'violet', 'gold']
        labels = list(counts.keys())
        values = list(counts.values())

        #Pie chart skips bases with zero count so their labels don't overlap
        pie_labels = []
        pie_values = []
        pie_colors = []
        for i in range(4):
            if values[i] > 0:
                pie_labels.append(labels[i])
                pie_values.append(values[i])
                pie_colors.append(colors[i])

        col1, col2 = st.columns(2)

        with col1:
            fig1, ax1 = plt.subplots()
            ax1.pie(pie_values, labels=pie_labels, colors=pie_colors, startangle=90, autopct='%1.1f%%')
            ax1.set_title("Nucleotide Composition - Pie Chart")
            ax1.axis('equal')
            st.pyplot(fig1)
            plt.close(fig1)

        with col2:
            fig2, ax2 = plt.subplots()
            ax2.bar(labels, values, color=colors)
            ax2.set_title("Nucleotide Counts - Bar Chart")
            ax2.set_xlabel("Nucleotide")
            ax2.set_ylabel("Count")
            ax2.grid(axis='y')
            st.pyplot(fig2)
            plt.close(fig2)
else:
    st.info("No sequence found yet. Paste a sequence or upload a FASTA file to begin.")