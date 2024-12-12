from flask import Flask, render_template, request, session, redirect, url_for
from Bio import SeqIO
import os
import sys
sys.path.append(os.path.abspath('../programa'))
from aliniacion_seqüencias import need_wunch_MM, need_wunch_MT, smith_waterman_MM, smith_waterman_MT
from Matrices import blosum62, pam250

app = Flask(__name__)
app.secret_key = 'asd'

# Folder to save uploaded files
UPLOAD_FOLDER = 'uploads'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# Allowed extensions
ALLOWED_EXTENSIONS = {'fasta'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

#pag principal
@app.route('/')
def home():
    return render_template('index.html')

# si aliniamiento global, se sube fasta y se decide si es DNA o Prot
@app.route('/Global', methods=['POST', 'GET'])
def upload_file():
    if request.method == 'POST':
        if 'file' not in request.files:
            return 'No file part in the request', 400
        file = request.files['file']
        if file.filename == '':
            return 'No file selected', 400
        if file and allowed_file(file.filename):
            filename = file.filename
            session['uploaded_file'] = filename
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

            molecula = request.form.get('opciones')
            if molecula == 'dna':
                return redirect(url_for('choose_parameters'))
            elif molecula == 'prot':
                return redirect(url_for('choose_parameters2'))
    return render_template('Global.html')

# si aliniamiento local, se sube fasta y se dice si es DNA o Prot
@app.route('/Local', methods=['POST', 'GET'])
def upload_file_1():
    if request.method == 'POST':
        if 'file' not in request.files:
            return 'No file part in the request', 400
        file = request.files['file']
        if file.filename == '':
            return 'No file selected', 400
        if file and allowed_file(file.filename):
            filename = file.filename
            session['uploaded_file'] = filename
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)

            molecula = request.form.get('opciones')
            if molecula == 'dna':
                return redirect(url_for('choose_parameters1'))
            elif molecula == 'prot':
                return redirect(url_for('choose_parameters3'))
    return render_template('Local.html')

# escogemos parametros y llamamos a la funcion que alinea
@app.route('/GlobalDNA', methods=['POST', 'GET'])
def choose_parameters():
    if request.method == 'POST':
        filename = session.get('uploaded_file')
        if not filename:
            return 'No file found in session', 400

        ruta = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        if not os.path.exists(ruta):
            return f'File {filename} not found in uploads folder', 400

        seq = []
        id = []

        for sec in SeqIO.parse(ruta, "fasta"):
            seq.append(sec.seq)
            id.append(sec.id)
        #inputs de la funcion
        seq_one = seq[0]
        seq_two = seq[1]
        match = int(request.form.get('match'))
        mismatch = int(request.form.get('mismatch'))
        gap = int(request.form.get('gap'))
        #funcion alin
        seq1, seq2, score = need_wunch_MM(seq_one, seq_two, match, mismatch, gap)
        return render_template('Resultado.html', result=[seq1, seq2, score])
    return render_template('GlobalDNA.html')

@app.route('/LocalDNA', methods=['POST', 'GET'])
def choose_parameters1():
    if request.method == 'POST':
        filename = session.get('uploaded_file')
        if not filename:
            return 'No file found in session', 400

        ruta = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        if not os.path.exists(ruta):
            return f'File {filename} not found in uploads folder', 400

        seq = []
        id = []

        for sec in SeqIO.parse(ruta, "fasta"):
            seq.append(sec.seq)
            id.append(sec.id)
        #inputs de la funcion
        seq_one = seq[0]
        seq_two = seq[1]
        match = int(request.form.get('match'))
        mismatch = int(request.form.get('miss'))
        gap = int(request.form.get('gap'))
        num = int(request.form.get('num'))
        #funcion alin
        seq1, seq2, score = smith_waterman_MM(seq_one, seq_two, match, mismatch, gap, num)
        return render_template('Resultado.html', result=[seq1, seq2, score])
    return render_template('LocalDNA.html')

@app.route('/GlobalProt', methods=['POST', 'GET'])
def choose_parameters2():
    if request.method == 'POST':
        filename = session.get('uploaded_file')
        if not filename:
            return 'No file found in session', 400

        ruta = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        if not os.path.exists(ruta):
            return f'File {filename} not found in uploads folder', 400

        seq = []
        id = []

        for sec in SeqIO.parse(ruta, "fasta"):
            seq.append(sec.seq)
            id.append(sec.id)
        #inputs de la funcion
        seq_one = seq[0]
        seq_two = seq[1]

        if request.form.get('tipo') == 'manual':
            match = int(request.form.get('match'))
            mismatch = int(request.form.get('miss'))
            gap = int(request.form.get('gap'))
            # funcion alin
            seq1, seq2, score = need_wunch_MM(seq_one, seq_two, match, mismatch, gap)
            return render_template('Resultado.html', result=[seq1, seq2, score])
        else:
            if request.form.get('matrix') == 'Blosum':
                matriz = blosum62
            else:
                matriz = pam250
            gaps = int(request.form.get('gaps'))
            #funcion ali
            seq1, seq2, score = need_wunch_MT(seq_one, seq_two, matriz, gaps)
            return render_template('Resultado.html', result=[seq1, seq2, score])
    return render_template('GlobalProt.html')

@app.route('/LocalProt', methods=['POST', 'GET'])
def choose_parameters3():
    if request.method == 'POST':
        filename = session.get('uploaded_file')
        if not filename:
            return 'No file found in session', 400

        ruta = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        if not os.path.exists(ruta):
            return f'File {filename} not found in uploads folder', 400

        seq = []
        id = []

        for sec in SeqIO.parse(ruta, "fasta"):
            seq.append(sec.seq)
            id.append(sec.id)
        #inputs de la funcion
        seq_one = seq[0]
        seq_two = seq[1]

        if request.form.get('tipo') == 'manual':
            match = int(request.form.get('match'))
            mismatch = int(request.form.get('miss'))
            gap = int(request.form.get('gap'))
            num = int(request.form.get('num_out'))
            # funcion alin
            seq1, seq2, score = smith_waterman_MM(seq_one, seq_two, match, mismatch, gap, num)
            return render_template('Resultado.html', result=[seq1, seq2, score])
        else:
            if request.form.get('matrix') == 'Blosum':
                matriz = blosum62
            else:
                matriz = pam250
            gaps = int(request.form.get('gaps'))
            num = int(request.form.get('num_out'))
            #funcion ali
            seq1, seq2, score = smith_waterman_MT(seq_one, seq_two, matriz, gaps, num)
            return render_template('Resultado.html', result=[seq1, seq2, score])
    return render_template('LocalProt.html')

if __name__ == '__main__':
    app.run(debug = True)