var x; // Array, microstate sequence
var S0; // parsed set of states
var S1; // mapped set of states
var ns; // integer, number of (unique) symbols
var lmap = {}; // dict, map input labels and mapped integers

function parseInput(){
    if (document.getElementById("input-option").value == 'paste') {
        var s = document.getElementById("paste-input").value;
    }
    else if (document.getElementById("input-option").value == 'upload') {
        var s = document.getElementById('input_data').value;
    } 
    //var raw = s.replace(/[^\d-]/g, '').trim(); // integers only
    var raw = s.replace(/[^a-zA-z0-9]/g, '').trim(); // letters and integers
    //S0 = [...new Set(raw)] // unique values of x
    S0 = Array.from(new Set(raw)); // unique values of x
    S1 = S0.map(a => S0.indexOf(a)); // map to 0,1,...
    S0.forEach(s => lmap[S0.indexOf(s)] = s); //lmap['mapped idx']='parsed idx'
    ns = S1.length;
    x = Array.prototype.map.call(raw, a => S0.indexOf(a));
    let raw_arr = Array.prototype.map.call(raw, a => a);
    
    document.getElementById("output-sequence-in").innerHTML = 
    "<strong>Parsed sequence: </strong>"+raw_arr.slice(0,30).join(', ')+"...";
    
    document.getElementById("output-symbols-parsed").innerHTML = 
    "<strong>Parsed symbols: </strong>" + S0.join(', ');
    
    document.getElementById("output-symbols-mapped").innerHTML = 
    "<strong>Mapped symbols: </strong>" + S1.join(', ');
    
    document.getElementById("output-sequence-mapped").innerHTML = 
    "<strong>Mapped sequence: </strong>" + x.slice(0,30).join(', ') + "...";
    
    //plotData(x)
}

function log2(x) {
    return Math.log(x)/Math.log(2.0);
}

function analyze(){
    // symbol distribution
    var p_hat = getP(x,ns)
    // transition matrix
    var T_hat = getT(x,ns)
    drawTable(T_hat);
    // entropy
    var h0 = H_1(p_hat);
    var hmax = log2(ns);
    // entropy rate
    var kmax = 7;
    var h1 = entropyRate(x,ns,kmax);
    // active information storage (AIS)
    var h2 = ais(x,ns,kmax);
    // autoinformation function (AIF)
    var lmax = 100;
    //var mi = aif(x, ns, lmax);
    //plotAif(aif);
    var mip = paif(x,ns,kmax);
    
    //var y = surrogate_mc(p_hat, T_hat, ns, x.length);
    
    var alpha = 0.05;
    
    var pval0 = testMarkov0(x, ns, alpha);
    var txtMarkov0 = "difference between Markov-0 and Markov-1" + 
    ((pval0 < alpha) ? " is " : " is not ") + "significant" + 
    " (p=" + pval0.toFixed(4) + ")";
    
    var pval1 = testMarkov1(x, ns, alpha);
    var txtMarkov1 = "difference between Markov-1 and Markov-2" + 
    ((pval1 < alpha) ? " is " : " is not ") + "significant" + 
    " (p=" + pval1.toFixed(4) + ")";
    
    // write output:
    console.log(lmap);
    document.getElementById("output-dist").innerHTML = 
    "<strong>Distribution: </strong>" + 
    "p(" + S0.join(', ') + ") = " + 
    p_hat.map((p,i) => p_hat[i].toFixed(2)).join(', ');
    /*
    "<strong>Distribution: </strong>" + 
    p_hat.map((val,idx) => "\np(" + S0[idx] + ") = " + val.toFixed(2));
    */
    
    document.getElementById("output-markov0").innerHTML = 
    "<strong>Markov-0 test: </strong>" + txtMarkov0;
    
    document.getElementById("output-markov1").innerHTML = 
    "<strong>Markov-1 test: </strong>" + txtMarkov1;
    
    document.getElementById("output-entropy").innerHTML = 
    "Result: " + h0.toFixed(3) + " bits" + 
    " (max. entropy: " + hmax.toFixed(3) + " bits)";
    
    document.getElementById("output-entropy-rate").innerHTML = 
    "Result for k = " + kmax + ": " + h1.toFixed(3) + " bits/sample";
    
    document.getElementById("output-ais").innerHTML = 
    "Result for k = 1..." + kmax + ": " +
    h2.map(a => a.toFixed(3)).join(', ');
    
    document.getElementById("output-paif").innerHTML = 
    "Result for k = 1..." + kmax + ": " +
    mip.map(a => a.toFixed(3)).join(', ');
}

function getP(y,n){
    var p = Array(n).fill(0);
    for (var i=0; i<y.length; i++) {p[y[i]]++}
    for (var j=0; j<p.length; j++) {p[j] /= y.length}
    return p;
}

function getT(y,n){
    var T = [];
    for (var i=0; i<n; i++) {T.push(Array(n).fill(0))}
    for (var t=0; t<y.length-1; t++) { T[y[t]][y[t+1]]++ }
    var rsum = Array(n).fill(0); // row sum
    for (var i=0; i<n; i++) {
        for (var j=0; j<n; j++) {
            rsum[i] += T[i][j];
        }
    }
    for (var i=0; i<n; i++) {
        if (rsum[i] > 0) {
            for (var j=0; j<n; j++) {
                T[i][j] /= rsum[i];
            }
        }
    }
    return T;
}

function H_1(p){
    var h = 0.0;
    for (var i=0; i<p.length; i++) {
        if (p [i] > 0) {
            h -= (p[i]*log2(p[i]));
        }
    }
    return h;
}

function H_k(x, ns, k) {
    // joint entropy, dim k
    let nmax = x.length-k;
    let p = {}; // joint distribution as dict, keys=data tuples
    for (var t=0; t<nmax; t++) {
        let y = x.slice(t,t+k);
        if (y in p) { p[y]++; }
        else { p[y] = 1.0; }
    }
    for (y in p) { p[y] /= nmax; }
    let hk = 0.0;
    for (y in p) { hk -= (p[y]*log2(p[y])); }
    return hk
}

function entropyRate(x,ns,kmax) {
    var h = Array(kmax).fill(0);
    var ks = Array(kmax);
    for (var k=0; k<kmax; k++) {
        h[k] = H_k(x, ns, k+1) // joint entropy for history length k+1
        ks[k] = k+1; // history lengths
    }
    // polynomial fit
    var lfit = linfit(h,ks);
    var a = lfit['a'];
    var b = lfit['b'];
    /*
    console.log("ks: " + ks);
    console.log("h: " + h);
    console.log("a: " + a);
    console.log("b: " + b);
    */
    //plotEntropyRate(ks,h,a,b)
    return a
}

function plotEntropyRate(ks,h,a,b) {
    var fit = Array(ks.length) //.fill(0);a*ks+b
    for (var i=0; i<ks.length; i++) { fit[i] = a*ks[i] + b; }
    var trace1 = {
        x: ks,
        y: h,
        type: 'scatter'
    };
    var trace2 = {
        x: ks,
        y: fit,
        type: 'line'
    };
    var layout = {
        xaxis: {autorange: true},
        yaxis: {autorange: true}
    };
    var data = [trace1, trace2];
    Plotly.newPlot('plot', data, layout);
}

function linfit(x,y){
        var fit = {};
        var n = x.length;
        var sum_x = 0;
        var sum_y = 0;
        var sum_xy = 0;
        var sum_xx = 0;
        var sum_yy = 0;
        for (var i=0; i<x.length; i++) {
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += (x[i]*y[i]);
            sum_xx += (x[i]*x[i]);
            sum_yy += (y[i]*y[i]);
        }
        // y_ = a*x + b
        fit['a'] = (n*sum_xy - sum_x*sum_y)/(n*sum_yy - sum_y*sum_y);
        fit['b'] = (sum_x - fit.a*sum_y)/n;
        //lr['r2'] = Math.pow((n*sum_xy - sum_x*sum_y)/Math.sqrt((n*sum_xx-sum_x*sum_x)*(n*sum_yy-sum_y*sum_y)),2);
        return fit
}

function aif(x, ns, lmax) {
    let n = x.length
    var mi = Array(lmax).fill(0);
    for (var l=0; l<lmax; l++) {
        let nmax = n - l;
        let p1 = Array(ns).fill(0);
        let p2 = Array(ns).fill(0);
        let p12 = [];
        for (var i=0; i<n; i++) { p12.push(Array(ns).fill(0)); }
        for (i=0; i<nmax; i++) {
            let i1 = x[i];
            let i2 = x[i+l];
            p1[i1]++;
            p2[i2]++;
            p12[i1][i2]++;
        }
        for (var i=0; i<ns; i++) {
            p1[i] /= nmax;
            p2[i] /= nmax;
            for (var j=0; j<ns; j++) {
                p12[i][j] /= nmax;
            }
        }
        // compute entropies
        let H1 = 0.0; 
        let H2 = 0.0;
        let H12 = 0.0;
        for (var i=0; i<ns; i++) {
            if (p1[i] > 0) { H1 -= (p1[i]*log2(p1[i])); }
            if (p2[i] > 0) { H2 -= (p2[i]*log2(p2[i])); }
            for (var j=0; j<ns; j++) {
                if (p12[i][j] > 0) { H12 -= (p12[i][j]*log2(p12[i][j])); }
            }
        }
        mi[l] = H1 + H2 - H12;
    }
    return mi
}

function paif(x, ns, kmax) {
    let n = x.length
    var y = Array(kmax).fill(0);
    var aif_ = aif(x,ns,2);
    y[0] = aif_[0];
    y[1] = aif_[1];
    var h1, h2, h3;
    for (k=2; k<kmax; k++) {
        h1 = H_k(x,ns,k);
        h2 = H_k(x,ns,k-1);
        h3 = H_k(x,ns,k+1);
        y[k] = 2*h1 - h2 - h3;
    }
    return y
}

function ais(x, ns, kmax) {
    let n = x.length;
    var y = Array(kmax).fill(0);
    for (var k=0; k<kmax; k++) {
        y[k] = H_k(x,ns,1) + H_k(x,ns,k) - H_k(x,ns,k+1);
    }
    return y
}

function surrogate_mc(p, T, ns, n) {
    /* Surrogate MC using the eq. distribution p, transition matrix T,
    ns symbols and length n */
    let r = Math.random(); // ~U[0,1], random threshold
    let s = 0.0;
    let y = p[s];
    while ( y < r ) {
        s++;
        y += p[s];
    }
    let Y = new Array();
    Y.push(s); // initial state according to p
    // iterate ...
    for (var i=1; i<n; i++) {
        r = Math.random(); // ~U[0,1], random threshold
        s = Y[i-1]; // currrent state
        let t = 0; // trial state
        y = T[s][t]; // transition rate to trial state
        while ( y < r ) {
            t++; // next trial state
            y += T[s][t];
        }
        Y.push(t); // store current state
    }
    return Y
}

function testMarkov0(x, ns, alpha) {
    /* Test zero-order Markovianity of symbolic sequence x with ns symbols
    Null hypothesis: zero-order MC (iid) <=>
    p( X[t] ) , p( X[t+1] ) independent, cf. Kullback, Technometrics (1962)
    */
    n = x.length;
    let f_i = Array(ns).fill(0);
    let f_j = Array(ns).fill(0);
    let f_ij = [];
    for (var i=0; i<ns; i++) { f_ij.push(Array(ns).fill(0)); }
    // calculate f_ij p( X[t]=i, p( X[t+1]=j ) )
    for (var t=0; t<(n-1); t++) {
        let i = x[t];
        let j = x[t+1];
        f_ij[i][j]++;
        f_i[i]++;
        f_j[j]++;
    }
    let T = 0.0 // statistic
    for (var i=0; i<ns; i++) {
        for (var j=0; j<ns; j++) {
            let f = f_ij[i][j]*f_i[i]*f_j[j];
            if ( f != 0.0 ) {
                T += (f_ij[i][j]*Math.log((n*f_ij[i][j])/(f_i[i]*f_j[j])));
            }
        }
    }
    T *= 2.0;
    let df = (ns-1.0) * (ns-1.0);
    let p = 1.0 - Gammacdf(T/2,df/2);
    //console.log(T + ", " + df + ", " + p)
    return p
}

function testMarkov1(x, ns, alpha) {
    /* Test first-order Markovianity of symbolic sequence X with ns symbols
    Null hypothesis: first-order MC <=>
    p( X[t+1] | X[t] ) = p( X[t+1] | X[t], X[t-1] ),
    cf. Kullback, Technometrics (1962), Tables 8.1, 8.2, 8.6.'''
    */
    let n = x.length;
    let f_j = Array(ns).fill(0);
    let f_ij = [];
    for (var i=0; i<ns; i++) { f_ij.push(Array(ns).fill(0)); }
    let f_jk = [];
    for (var i=0; i<ns; i++) { f_jk.push(Array(ns).fill(0)); }
    let f_ijk = [];
    for (var i=0; i<ns; i++) { f_ijk.push([]); }
    for (var i=0; i<ns; i++) {
        for (var j=0; j<ns; j++) { 
            f_ijk[i].push(Array(ns).fill(0));
        }
    }
    for (var t=0; t<(n-2); t++) {
        let i = x[t];
        let j = x[t+1];
        let k = x[t+2];
        f_ijk[i][j][k]++;
        f_ij[i][j]++;
        f_jk[j][k]++;
        f_j[j]++;
    }
    let T = 0.0;
    for (var i=0; i<ns; i++) {
        for (var j=0; j<ns; j++) {
            for (var k=0; k<ns; k++) {
                f = f_ijk[i][j][k]*f_j[j]*f_ij[i][j]*f_jk[j][k];
                if ( f != 0.0 ) {
                    T += (f_ijk[i][j][k]*Math.log((f_ijk[i][j][k]*f_j[j])/(f_ij[i][j]*f_jk[j][k])))
                }
            }
        }
    }
    T *= 2.0;
    let df = ns*(ns-1)*(ns-1);
    let p = 1.0 - Gammacdf(T/2,df/2);
    // console.log(T + ", " + df + ", " + p)
    return p
}

function drawTable(T) {
    //document.getElementById('div1').innerHTML = "";
    nr = T.length;
    nc = T[0].length;
    var div1 = document.getElementById('div1');
    var tbl = document.createElement("table");
    for (var r=0; r<nr; r++) {
        var row = document.createElement("tr");
        for (var c=0; c<nc; c++) {
            var cell = document.createElement("td");
            var cellText = document.createTextNode(T[r][c].toFixed(2));
            cell.appendChild(cellText);
            row.appendChild(cell);
        }
        tbl.appendChild(row);
    }
    div1.appendChild(tbl);
}

function plotData(z){
    var trace1 = {
        x: [...Array(z.length).keys()],
        y: z,
        type: 'scatter'
    };
    var data = [trace1];
    Plotly.newPlot('plot', data);
}

function plotAif(z){
    var trace1 = {
        x: [...Array(z.length).keys()],
        y: z,
        type: 'scatter'
    };
    var layout = {
        xaxis: {type: 'lin', autorange: true},
        yaxis: {type: 'log', autorange: true}
    };
    var data = [trace1];
    Plotly.newPlot('plot', data, layout);
}

// adapted from: https://www.math.ucla.edu/~tom/distributions/chisq.html

function LogGamma(Z) {
    var c1 = 76.18009173;
    var c2 = -86.50532033;
    var c3 = 24.01409822;
    var c4 = -1.231739516;
    var c5 = 0.00120858003;
    var c6 = -0.00000536382;
    var S = 1 + c1/Z + c2/(Z+1) + c3/(Z+2) + c4/(Z+3) + c5/(Z+4) + c6/(Z+5);
    var c7 = 2.50662827465;
    var LG = (Z-0.5)*Math.log(Z+4.5) - (Z+4.5) + Math.log(S*c7);
	return LG
}

function Gcf(X,A) {
    // Good for X>A+1
    var A0 = 0;
	var B0 = 1;
	var A1 = 1;
	var B1 = X;
	var AOLD = 0;
	var N = 0;
	while (Math.abs((A1-AOLD)/A1) > 0.00001) {
		AOLD = A1;
		N = N+1;
		A0 = A1 + (N-A)*A0;
		B0 = B1 + (N-A)*B0;
		A1 = X*A0 + N*A1;
		B1 = X*B0 + N*B1;
		A0 = A0/B1;
		B0 = B0/B1;
		A1 = A1/B1;
		B1 = 1;
	}
	var Prob = Math.exp(A*Math.log(X)-X-LogGamma(A))*A1;
	return 1-Prob
}

function Gser(X,A) {
    // Good for X<A+1.
	var T9 = 1/A;
	var G = T9;
	var I = 1;
	while (T9 > G*0.00001) {
		T9 = T9*X/(A+I);
		G = G+T9;
		I = I+1;
	}
	G=G*Math.exp(A*Math.log(X)-X-LogGamma(A));
    return G
}

function Gammacdf(x,a) {
	var GI;
	if (x <= 0) {
		GI = 0;
	} else if (x < a+1) {
		GI = Gser(x,a);
	} else {
		GI = Gcf(x,a);
	}
	return GI
}
