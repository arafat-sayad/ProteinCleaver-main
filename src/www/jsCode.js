shinyjs.createmolart = function(unipID, seque, startPos, endPos, identifiable){
  
  document.getElementById("molart-viewer").remove();
  var divEl = document.createElement("div");
  divEl.setAttribute("id", "molart-viewer");
  document.getElementById("fluidRow-molart-viewer").appendChild(divEl);
  
  
  document.getElementById("molart-viewer").innerHTML = "";
  document.getElementById("molart-viewer").style.height="auto";
  
  var startPosArr = startPos.split(",");
  var endPosArr = endPos.split(",");
  var identifiableArr = identifiable.split(",");
  var feat = "";
  
  for (var i=0; i < startPosArr.length; i++) {
    if (identifiableArr[i] === '1') {
      feat = feat.concat("{type: 'IDENTIFIABLE PEPTIDE', category: 'IDENTIFIABLE PEPTIDES', description: 'Identifiable peptide', begin: ",startPosArr[i],", end: ",endPosArr[i],", color: '#22478a'},");
    }
    else if (identifiableArr[i] === '0') {
      feat = feat.concat("{type: 'UNDETECTABLE PEPTIDE', category: 'UNDETECTABLE PEPTIDES', description: 'Non identifiable peptide', begin: ",startPosArr[i],", end: ",endPosArr[i],", color: '#d12741'},");
    }
  }
  
  feat = feat.substring(0, feat.length - 1);
  feat = ("[").concat(feat);
  feat = feat.concat("]");
  var molArtFeatures = eval("(" + feat + ")");
  //console.log(molArtFeatures);
  
  molart = new MolArt({
    uniprotId: unipID,
    containerId: "molart-viewer",
    customDataSources: [
        {
            source: 'RANDOM',
            useExtension: false,
            data: {
                sequence: seque,
                features: molArtFeatures
            }
        }
    ]
  });
    
  var x = document.getElementById("showStructBtn");
  x.onclick = function() { shinyjs.togglemolart(); };

  var structBtn = document.getElementById("showStructBtn");
	structBtn.setAttribute("class", "btn btn-danger btn-sm");
	structBtn.innerHTML = "Hide protein features";
    
};

shinyjs.togglemolart = function() {
  var x = document.getElementById("molart-viewer");
  if (x.style.display === "none") {
    x.style.display = "block";
    var structBtn = document.getElementById("showStructBtn");
	  structBtn.setAttribute("class", "btn btn-danger btn-sm");
	  structBtn.innerHTML = "Hide protein features";
  } else {
    x.style.display = "none";
    var structBtn2 = document.getElementById("showStructBtn");
	  structBtn2.setAttribute("class", "btn btn-custom btn-sm");
	  structBtn2.innerHTML = "Show protein features";
  }
};


shinyjs.hidesequeNmolart = function() {
  document.getElementById("sequence-viewer").remove();
  var divEl1 = document.createElement("div");
  divEl1.setAttribute("id", "sequence-viewer");
  document.getElementById("fluidRow-sequence-viewer").appendChild(divEl1);
  
  document.getElementById("molart-viewer").remove();
  var divEl2 = document.createElement("div");
  divEl2.setAttribute("id", "molart-viewer");
  document.getElementById("fluidRow-molart-viewer").appendChild(divEl2);
};


shinyjs.seque = function(unipID, seque, geneName, protName, startPos, endPos, identifiable, unq) {
  
  document.getElementById("molart-viewer").innerHTML = "";
  document.getElementById("molart-viewer").style.height="auto";

  var seq=new Sequence(seque);
  var concat_title = unipID.concat(" | ", geneName.concat(" | ", protName));
  seq.render("#sequence-viewer", {"charsPerLine": 100, "toolbar": false, "title": concat_title, "search" : true}   );
  
  /*var onclickFun = function(e) {

  }; */

  var legend = [
    {name: "Identifiable peptide", color: "#22478a", underscore: false},
    {name: "Undetectable peptide", color: "#d12741", underscore: false},
    {name: "Cleavage site", color: "#CCCC99", underscore: false},
    {name: "Unique peptide", underscore: true}
    ];
  seq.addLegend(legend);
  
  
  var startPosArr = startPos.split(",");
  var endPosArr = endPos.split(",");
  var identifiableArr = identifiable.split(",");
  var identifiableCol = identifiable.replace(/1/gi, "#22478a");
  identifiableCol = identifiableCol.replace(/0/gi, "#d12741");
  var identifiableColArr = identifiableCol.split(",");
  var uniqArr = unq.split(",");
  var unqUnderSc = unq.replace(/1/gi, true);
  unqUnderSc = unqUnderSc.replace(/0/gi, false);
  var unqUnderScArr = unqUnderSc.split(",");
  var highlightColor = "#CCCC99";
  var seqCoverage = "";
  
  //for peptides (start to stop-1)
  for (var i=0; i < startPosArr.length; i++) {
    seqCoverage = seqCoverage.concat("{start:",startPosArr[i]-1,",end:",endPosArr[i]-1,",color:",JSON.stringify(identifiableColArr[i]),", underscore: ",unqUnderScArr[i],"}," );
  }
  
  //highlight cleavage sites
 for (var j=0; j < endPosArr.length; j++) {
   if (j !=  endPosArr.length - 1) { //highlight cleavage sites
     seqCoverage = seqCoverage.concat("{start:",endPosArr[j]-1,",end:",endPosArr[j],",color:",JSON.stringify(identifiableColArr[j]),",underscore: ",unqUnderScArr[j],", bgcolor: highlightColor}," );
   }
   else { //remove highlightColor in the last amino acid
     seqCoverage = seqCoverage.concat("{start:",endPosArr[j]-1,",end:",endPosArr[j],",color:",JSON.stringify(identifiableColArr[j]),",underscore: ",unqUnderScArr[j],"}," );
   }
  } 

  seqCoverage = seqCoverage.substring(0, seqCoverage.length - 1);
  seqCoverage = ("[").concat(seqCoverage);
  seqCoverage = seqCoverage.concat("]");
  
  var sequenceCoverage = eval("(" + seqCoverage + ")");
  seq.coverage(sequenceCoverage);  
  
  var aEl = document.createElement("a");
  aEl.setAttribute("id", "showStructBtn");
	aEl.href = "#molart-viewer";
	aEl.onclick = function() { shinyjs.createmolart(unipID, seque, startPos, endPos, identifiable); };
	aEl.className = "btn btn-custom btn-sm";
	aEl.setAttribute("style", "float: right; vertical-align:top;");
	
	var spanEl = document.createElement("span");
  spanEl.classList.add("spanVizBtn");
  aEl.appendChild(spanEl);
  aEl.insertAdjacentText("beforeend","Show protein features");
  
  var divEl = document.getElementsByClassName("coverageLegend")[0];
  divEl.appendChild(aEl);
  
  /* seq.onMouseSelection(function(elem){
        console.log(elem.detail);
    }
  );
  seq.onSubpartSelected(function(elem){
        console.log(elem.detail);
    }
  ); */

};


shinyjs.rezet = function() {
  //history.go(0);
  location.reload(true);
};


shinyjs.swalErrorAlert = function() {
    Swal.fire({
        title: "Something went wrong",
        //html:  missingColumn,
        html: "An error occured while uploading your data. Please check your input file and try again",
        icon: "error",
        width: 700,
        padding: '2em',
        showConfirmButton: true,
        allowEnterKey: true,
        allowEscapeKey: true,
        showCloseButton: true,
        allowOutsideClick: true
    });
};


shinyjs.swalAlert = function(header, message, iconType) {
    Swal.fire({
        title: header,
        html: message,
        icon: iconType,
        width: 700,
        padding: '2em',
        showConfirmButton: false,
        allowEnterKey: true,
        allowEscapeKey: true,
        showCloseButton: true,
        allowOutsideClick: true
    });
};


shinyjs.disableTab = function(name) {
var tab = $(".nav li a[data-value=" + name + "]");
  tab.bind("click.tab", function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass("disabled");
};


shinyjs.enableTab = function(name) {
  var tab = $(".nav li a[data-value=" + name + "]");
  tab.unbind("click.tab");
  tab.removeClass("disabled");
};


$(document).on('click', 'tr', function( e, dt, type, cell, originalEvent ){
    $(this).addClass('selected');
});
