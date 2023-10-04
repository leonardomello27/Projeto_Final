require(["dojo/ready"], function(ready){
  ready(function(){
    dojo.style(dojo.byId('slwebview_InspectorPane'), "display", "none");
    dojo.style(dojo.byId('slwebview_InspectorPane_splitter'), "display", "none");
    dojo.style(dojo.byId('dijit__TemplatedMixin_1'), "display", "none");
    dijit.byId('slwebview').resize();
  });
});